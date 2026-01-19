import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import random

# --- 1. Define the Super-Resolution Model ---
# A simple 3D CNN that takes a coarse SDF and outputs a fine SDF.
# We use a residual learning approach: Output = Input_Upscaled + Residual
class SDFSuperResNet(nn.Module):
    def __init__(self, scale_factor=2):
        super(SDFSuperResNet, self).__init__()
        self.scale_factor = scale_factor
        
        # Upsampling layer (Trilinear interpolation)
        self.upsample = nn.Upsample(scale_factor=scale_factor, mode='trilinear', align_corners=True)
        
        # Feature extraction (Residual blocks)
        self.conv1 = nn.Conv3d(1, 64, kernel_size=3, padding=1)
        self.relu = nn.ReLU(inplace=True)
        
        self.res_blocks = nn.Sequential(
            nn.Conv3d(64, 64, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 64, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 64, kernel_size=3, padding=1),
            nn.ReLU(inplace=True)
        )
        
        self.conv_final = nn.Conv3d(64, 1, kernel_size=3, padding=1)

    def forward(self, x):
        # x shape: [Batch, 1, D, H, W] (Coarse grid)
        
        # 1. Base prediction: Simple interpolation
        x_up = self.upsample(x)
        
        # 2. Learn the residual (the details)
        out = self.relu(self.conv1(x_up))
        out = self.res_blocks(out)
        residual = self.conv_final(out)
        
        # 3. Combine
        return x_up + residual

# --- 2. Synthetic Data Generator ---
# Generates random SDF patches (Spheres and Boxes)
def generate_sdf_batch(batch_size, size=32):
    # Returns [Batch, 1, size, size, size]
    data = torch.zeros(batch_size, 1, size, size, size)
    
    grid = torch.stack(torch.meshgrid(
        torch.arange(size), torch.arange(size), torch.arange(size), indexing='ij'
    )) # [3, size, size, size]
    
    for i in range(batch_size):
        # Randomly choose primitive
        shape_type = random.choice(['sphere', 'box'])
        center = torch.rand(3) * size
        
        if shape_type == 'sphere':
            radius = random.uniform(size * 0.2, size * 0.45)
            dist = torch.sqrt(((grid - center.view(3, 1, 1, 1))**2).sum(0))
            sdf = dist - radius
            
        elif shape_type == 'box':
            dims = torch.rand(3) * (size * 0.4) + (size * 0.1)
            # Box SDF: length(max(d,0)) + min(max(d.x,max(d.y,d.z)),0.0)
            # where d = abs(p) - b
            p = grid - center.view(3, 1, 1, 1)
            d = torch.abs(p) - dims.view(3, 1, 1, 1)
            inside = torch.min(torch.max(d[0], torch.max(d[1], d[2])), torch.tensor(0.0))
            outside = torch.sqrt((torch.max(d, torch.tensor(0.0))**2).sum(0))
            sdf = inside + outside
            
        data[i, 0] = sdf
        
    return data

# --- 3. Training Loop ---
def train():
    # Hyperparameters
    SCALE = 2
    FINE_SIZE = 32
    COARSE_SIZE = FINE_SIZE // SCALE
    BATCH_SIZE = 16
    ITERATIONS = 1000
    
    # Set device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    model = SDFSuperResNet(scale_factor=SCALE).to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-4)
    criterion = nn.MSELoss()
    
    model.train()
    
    for it in range(ITERATIONS):
        # 1. Generate High-Res Ground Truth
        fine_sdf = generate_sdf_batch(BATCH_SIZE, FINE_SIZE).to(device)
        
        # 2. Generate Low-Res Input (Downsample)
        # ViennaLS lsResample picks the point with the smallest absolute value.
        # We simulate this using MaxPool3d on -abs(sdf) with indices to retrieve the original value (with sign).
        neg_abs_fine = -torch.abs(fine_sdf)
        _, indices = F.max_pool3d(neg_abs_fine, kernel_size=SCALE, stride=SCALE, return_indices=True)
        
        # Gather the original values using the indices
        fine_sdf_flat = fine_sdf.view(fine_sdf.size(0), fine_sdf.size(1), -1)
        indices_flat = indices.view(indices.size(0), indices.size(1), -1)
        coarse_sdf = torch.gather(fine_sdf_flat, 2, indices_flat).view(indices.size())
        
        coarse_sdf /= SCALE # Normalize SDF to coarse grid units
        
        # 3. Forward Pass
        optimizer.zero_grad()
        predicted_fine = model(coarse_sdf)
        
        # 4. Loss
        # We care most about the interface (SDF ~ 0). 
        # Target is fine_sdf normalized to coarse grid units (to match input scale)
        target = fine_sdf / SCALE
        weight = torch.exp(-torch.abs(target)) # Higher weight where SDF is close to 0
        loss = (weight * (predicted_fine - target)**2).mean()
        
        loss.backward()
        optimizer.step()
        
        if it % 100 == 0:
            print(f"Iter {it}: Loss = {loss.item():.6f}")
            
    print("Training finished.")
    
    # --- 4. Export Model ---
    # We export to TorchScript so it can be loaded in C++
    model.eval()
    example_input = torch.rand(1, 1, COARSE_SIZE, COARSE_SIZE, COARSE_SIZE).to(device)
    traced_script_module = torch.jit.trace(model, example_input)
    traced_script_module.save("sdf_super_res.pt")
    print("Model saved to sdf_super_res.pt")

if __name__ == "__main__":
    train()
# ```

# ### 3. How to use this in C++
# Once you run the python script and get `sdf_super_res.pt`, you can load it in your C++ application using LibTorch (the C++ frontend for PyTorch).

# In your `lsNeuralAdvect` callback:

# ```cpp
# // Inside your main C++ file
# #include <torch/script.h> // LibTorch header

# // Load model once
# torch::jit::script::Module module;
# try {
#    module = torch::jit::load("sdf_super_res.pt");
# } catch (const c10::Error& e) {
#    std::cerr << "Error loading model\n";
# }

# // Setup Neural Advection
# auto nnAdvect = ls::NeuralAdvect<double, 3>(myDomain, myVelocities);
# nnAdvect.setCoarseningFactor(2.0);

# nnAdvect.setSuperResolutionCallback([&module](auto coarseLS, auto fineLS) {
#    // 1. Convert coarseLS grid to Torch Tensor
#    // (Iterate grid, copy values to tensor)
#    // Tensor input = ...; // Shape {1, 1, D, H, W}

#    // 2. Run Inference
#    // std::vector<torch::jit::IValue> inputs;
#    // inputs.push_back(input);
#    // at::Tensor output = module.forward(inputs).toTensor();

#    // 3. Write output Tensor back to fineLS
#    // (Iterate tensor, write values to fineLS grid)
# });
