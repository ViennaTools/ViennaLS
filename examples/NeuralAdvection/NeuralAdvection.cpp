#include <iostream>
#include <vector>

// ViennaLS Includes
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsNeuralAdvect.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsCompareChamfer.hpp>
#include <vcTimer.hpp>

// LibTorch Include
#include <torch/script.h>
#include <torch/cuda.h>

using namespace viennals;

// Define a simple velocity field (e.g., growing)
class ConstantVelocity : public VelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> &, int,
                           const std::array<double, 3> &,
                           unsigned long) override {
    return 0.5; // Constant growth speed
  }
};

int main() {
  // 1. Setup Simulation Domain
  constexpr int D = 3;
  double gridDelta = 0.1;
  double bounds[6] = {-10, 10, -10, 10, -10, 10};
  Domain<double, D>::BoundaryType boundaryCons[D];
  for (int i = 0; i < D; ++i)
    boundaryCons[i] = Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet = Domain<double, D>::New(bounds, boundaryCons, gridDelta);
  double minCorner[3] = {-2.0, -2.0, -2.0};
  double maxCorner[3] = {2.0, 2.0, 2.0};
  MakeGeometry<double, D>(levelSet, Box<double, D>::New(minCorner, maxCorner)).apply();

  double advectionTime = 2.0;
  // Final box: Growth of 0.5 * 2.0 = 1.0 per side
  double finalMin[3] = {-3.0, -3.0, -3.0};
  double finalMax[3] = {3.0, 3.0, 3.0};

  // 2. Load the TorchScript Model
  torch::jit::script::Module model;
  try {
    // Ensure this path points to your trained model file
    model = torch::jit::load("sdf_super_res.pt");

    if (torch::cuda::is_available()) {
      std::cout << "LibTorch: CUDA available! Moving model to GPU."
                << std::endl;
      model.to(torch::kCUDA);
    } else {
      std::cout << "LibTorch: CUDA not available. Using CPU." << std::endl;
    }
  } catch (const c10::Error &e) {
    std::cerr << "Error loading the model: " << e.what() << std::endl;
    return -1;
  }

  double scaleFactor = 2.0;
  
  auto callback = [&](SmartPointer<Domain<double, D>> coarseLS,
                      SmartPointer<Domain<double, D>> fineLS) {
    std::cout << "  [Callback] Running Neural Super-Resolution..." << std::endl;

    auto &domain = coarseLS->getDomain();

    // A. Determine Bounding Box of Active Region
    // (In a real app, you would tile the domain into smaller chunks here)
    int minIdx[D], maxIdx[D];
    for (int i = 0; i < D; ++i) {
      minIdx[i] = std::numeric_limits<int>::max();
      maxIdx[i] = std::numeric_limits<int>::min();
    }

    viennahrle::ConstSparseIterator<Domain<double, D>::DomainType> it(domain);
    for (; !it.isFinished(); ++it) {
      if (it.isDefined()) {
        for (int i = 0; i < D; ++i) {
          if (it.getStartIndices(i) < minIdx[i]) minIdx[i] = it.getStartIndices(i);
          if (it.getStartIndices(i) > maxIdx[i]) maxIdx[i] = it.getStartIndices(i);
        }
      }
    }

    // Add padding for context
    int padding = 2;
    for (int i = 0; i < D; ++i) {
      minIdx[i] -= padding;
      maxIdx[i] += padding;
    }

    // Check if domain is empty or invalid
    bool isEmpty = false;
    for (int i = 0; i < D; ++i) {
      if (maxIdx[i] < minIdx[i]) isEmpty = true;
    }
    if (isEmpty) {
      std::cout << "  [Callback] Coarse domain is empty. Skipping." << std::endl;
      return;
    }

    int dims[D];
    for (int i = 0; i < D; ++i) dims[i] = maxIdx[i] - minIdx[i] + 1;

    std::cout << "  [Callback] Coarse Grid Bounds: [" << minIdx[0] << ", " << minIdx[1] << ", " << minIdx[2] << "] to ["
              << maxIdx[0] << ", " << maxIdx[1] << ", " << maxIdx[2] << "]" << std::endl;
    std::cout << "  [Callback] Tensor Dims: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;

    // B. Create Input Tensor (Batch=1, Channel=1, D, H, W)
    // Initialize with a safe background value (positive = outside)
    std::vector<int64_t> tensorShape = {1, 1};
    if constexpr (D == 3) {
      tensorShape.push_back(dims[2]);
      tensorShape.push_back(dims[1]);
      tensorShape.push_back(dims[0]);
    } else {
      tensorShape.push_back(1); // Depth 1 for 2D
      tensorShape.push_back(dims[1]);
      tensorShape.push_back(dims[0]);
    }
    torch::Tensor inputTensor = torch::full(tensorShape, 5.0, torch::kFloat32);
    auto accessor = inputTensor.accessor<float, 5>();

    int definedPointCount = 0;

    // Fill tensor from sparse grid
    viennahrle::ConstSparseIterator<Domain<double, D>::DomainType> it2(domain);
    for (; !it2.isFinished(); ++it2) {
      if (it2.isDefined()) {
        int x = it2.getStartIndices(0) - minIdx[0];
        int y = it2.getStartIndices(1) - minIdx[1];
        int z = 0;
        if constexpr (D == 3) z = it2.getStartIndices(2) - minIdx[2];

        bool inBounds = (x >= 0 && x < dims[0] && y >= 0 && y < dims[1]);
        if constexpr (D == 3) inBounds = inBounds && (z >= 0 && z < dims[2]);

        if (inBounds) {
          accessor[0][0][z][y][x] = static_cast<float>(it2.getValue());
          definedPointCount++;
        }
      }
    }
    std::cout << "  [Callback] Filled tensor with " << definedPointCount << " defined points." << std::endl;

    if (torch::cuda::is_available()) inputTensor = inputTensor.to(torch::kCUDA);

    // C. Run Inference
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(inputTensor);
    torch::Tensor outputTensor = model.forward(inputs).toTensor();
    outputTensor = outputTensor.to(torch::kCPU);

    // D. Write Output back to Fine Level Set
    auto outAccessor = outputTensor.accessor<float, 5>();
    int outDims[3] = {static_cast<int>(outputTensor.size(4)), 
                      static_cast<int>(outputTensor.size(3)), 
                      static_cast<int>(outputTensor.size(2))};

    std::vector<std::pair<viennahrle::Index<D>, double>> newPoints;
    
    // Calculate offset in fine grid coordinates
    int fineMinIdx[D];
    for (int i = 0; i < D; ++i) fineMinIdx[i] = minIdx[i] * static_cast<int>(scaleFactor);

    int insertedCount = 0;
    double minVal = 1e9, maxVal = -1e9;

    for (int z = 0; z < outDims[2]; ++z) {
      for (int y = 0; y < outDims[1]; ++y) {
        for (int x = 0; x < outDims[0]; ++x) {
          // The NN outputs distance in Coarse Grid Units.
          // We must multiply by scaleFactor to convert to Fine Grid Units.
          float val = outAccessor[0][0][z][y][x] * static_cast<float>(scaleFactor);

          if(val < minVal) minVal = val;
          if(val > maxVal) maxVal = val;

          // Only insert points near the interface to keep it sparse
          // Relaxed threshold to 4.0 to ensure we catch the interface
          if (std::abs(val) < 4.0) {
            viennahrle::Index<D> idx;
            idx[0] = fineMinIdx[0] + x;
            idx[1] = fineMinIdx[1] + y;
            if constexpr (D == 3) idx[2] = fineMinIdx[2] + z;
            newPoints.push_back({idx, static_cast<double>(val)});
            insertedCount++;
          }
        }
      }
    }

    std::cout << "  [Callback] NN Output Range: [" << minVal << ", " << maxVal << "]. Inserted " << insertedCount << " points." << std::endl;
    fineLS->insertPoints(newPoints);
  };

  // --- Generate Ideal Solution ---
  auto idealLS = Domain<double, D>::New(bounds, boundaryCons, gridDelta);
  MakeGeometry<double, D>(idealLS, Box<double, D>::New(finalMin, finalMax)).apply();
  
  // Generate Ideal Coarse Solution for fair comparison
  auto idealCoarseLS = Domain<double, D>::New(bounds, boundaryCons, gridDelta * scaleFactor);
  MakeGeometry<double, D>(idealCoarseLS, Box<double, D>::New(finalMin, finalMax)).apply();

  // --- 1. Fine Grid Simulation ---
  {
    std::cout << "\n--- Starting Fine Grid Simulation ---" << std::endl;
    auto fineLS = Domain<double, D>::New(bounds, boundaryCons, gridDelta);
    MakeGeometry<double, D>(fineLS, Box<double, D>::New(minCorner, maxCorner)).apply();
    
    auto velocities = SmartPointer<ConstantVelocity>::New();
    Advect<double, D> advect(fineLS, velocities);
    advect.setAdvectionTime(advectionTime);
    advect.setSpatialScheme(SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
    advect.setTemporalScheme(TemporalSchemeEnum::FORWARD_EULER);

    Timer timer;
    timer.start();
    advect.apply();
    timer.finish();
    std::cout << "Fine Simulation Time: " << timer.currentDuration / 1e9 << "s" << std::endl;

    CompareChamfer<double, D> chamfer(idealLS, fineLS);
    chamfer.apply();
    std::cout << "Fine Chamfer Distance: " << chamfer.getChamferDistance() << std::endl;

    auto mesh = SmartPointer<Mesh<double>>::New();
    ToSurfaceMesh<double, D>(fineLS, mesh).apply();
    VTKWriter<double>(mesh, "Result_Fine.vtp").apply();
  }

  // --- 2. Coarse Grid Simulation ---
  {
    std::cout << "\n--- Starting Coarse Grid Simulation ---" << std::endl;
    auto coarseLS = Domain<double, D>::New(bounds, boundaryCons, gridDelta * scaleFactor);
    MakeGeometry<double, D>(coarseLS, Box<double, D>::New(minCorner, maxCorner)).apply();
    
    auto velocities = SmartPointer<ConstantVelocity>::New();
    Advect<double, D> advect(coarseLS, velocities);
    advect.setAdvectionTime(advectionTime);
    advect.setSpatialScheme(SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
    advect.setTemporalScheme(TemporalSchemeEnum::FORWARD_EULER);

    Timer timer;
    timer.start();
    advect.apply();
    timer.finish();
    std::cout << "Coarse Simulation Time: " << timer.currentDuration / 1e9 << "s" << std::endl;

    CompareChamfer<double, D> chamfer(idealLS, coarseLS);
    chamfer.apply();
    std::cout << "Coarse Chamfer Distance: " << chamfer.getChamferDistance() << std::endl;

    auto mesh = SmartPointer<Mesh<double>>::New();
    ToSurfaceMesh<double, D>(coarseLS, mesh).apply();
    VTKWriter<double>(mesh, "Result_Coarse.vtp").apply();
  }

  // --- 3. Neural Advection Simulation ---
  {
    std::cout << "\n--- Starting Neural Advection Simulation ---" << std::endl;
    auto neuralLS = Domain<double, D>::New(bounds, boundaryCons, gridDelta);
    MakeGeometry<double, D>(neuralLS, Box<double, D>::New(minCorner, maxCorner)).apply();
    
    auto velocities = SmartPointer<ConstantVelocity>::New();
    NeuralAdvect<double, D> nnAdvect(neuralLS, velocities);
    nnAdvect.setCoarseningFactor(scaleFactor);
    nnAdvect.setSuperResolutionCallback(callback);

    Timer timer;
    timer.start();
    nnAdvect.apply(advectionTime);
    timer.finish();
    std::cout << "Neural Simulation Time: " << timer.currentDuration / 1e9 << "s" << std::endl;

    CompareChamfer<double, D> chamfer(idealLS, neuralLS);
    chamfer.apply();
    std::cout << "Neural Chamfer Distance: " << chamfer.getChamferDistance() << std::endl;

    auto mesh = SmartPointer<Mesh<double>>::New();
    ToSurfaceMesh<double, D>(neuralLS, mesh).apply();
    VTKWriter<double>(mesh, "Result_Neural.vtp").apply();
  }

  return 0;
}