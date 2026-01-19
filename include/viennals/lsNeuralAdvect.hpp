#pragma once

#include <lsAdvect.hpp>
#include <lsResample.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <functional>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>


namespace viennals {

/// Experimental class demonstrating Multi-Resolution Advection
/// accelerated by Neural Networks.
template <class T, int D> class NeuralAdvect {
  SmartPointer<Domain<T, D>> fineLevelSet;
  SmartPointer<VelocityField<T>> velocities;
  
  // Factor to downsample (e.g., 2.0 means grid becomes 2x coarser)
  double coarseningFactor = 2.0; 
  
  // Callback for the Neural Network Inference
  // Input: Coarse Domain, Output: Fine Domain (Super-Resolved)
  using SuperResolutionCallback = std::function<void(SmartPointer<Domain<T, D>>, SmartPointer<Domain<T, D>>)>;
  SuperResolutionCallback nnInference = nullptr;

public:
  NeuralAdvect(SmartPointer<Domain<T, D>> ls, SmartPointer<VelocityField<T>> vel) 
    : fineLevelSet(ls), velocities(vel) {}

  void setCoarseningFactor(double factor) { coarseningFactor = factor; }
  
  void setSuperResolutionCallback(SuperResolutionCallback cb) {
      nnInference = cb;
  }

  void apply(double advectionTime) {
      if (coarseningFactor <= 1.0) {
          // Fallback to standard advection if no coarsening
          Advect<T, D> standardAdvect(fineLevelSet, velocities);
          standardAdvect.setAdvectionTime(advectionTime);
          standardAdvect.apply();
          return;
      }

      // 1. DOWNSAMPLE
      // Create coarse domain
      double fineDelta = fineLevelSet->getGrid().getGridDelta();
      double coarseDelta = fineDelta * coarseningFactor;
      
      auto coarseLevelSet = SmartPointer<Domain<T, D>>::New(coarseDelta);
      
      // Use Resample tool
      Resample<T, D> downsampler(fineLevelSet, coarseLevelSet, coarseDelta);

      // Enable sub-grid correction 
      downsampler.setSubGridCorrection(true);
      downsampler.apply();

      // Ensure coarse level set has enough width for advection
      // This prevents the surface from moving out of the defined narrow band
      Expand<T, D>(coarseLevelSet, 6).apply();
      
      {
        auto mesh = Mesh<T>::New();
        ToMesh<T, D>(fineLevelSet, mesh).apply();
        VTKWriter<T>(mesh, "fineLevelSet.vtu").apply();
        ToSurfaceMesh<T, D>(fineLevelSet, mesh).apply();
        VTKWriter<T>(mesh, "fineLevelSet.vtp").apply();
      }
      {
        auto mesh = Mesh<T>::New();
        ToMesh<T, D>(coarseLevelSet, mesh).apply();
        VTKWriter<T>(mesh, "coarseLevelSet.vtu").apply();
        ToSurfaceMesh<T, D>(coarseLevelSet, mesh).apply();
        VTKWriter<T>(mesh, "coarseLevelSet.vtp").apply();
      }


      // 2. COARSE ADVECTION
      // Run standard ViennaLS advection on the coarse grid.
      // This is much faster because:
      // a) Fewer points (N / factor^D)
      // b) Larger time steps (dt_max ~ coarseDelta)
      Advect<T, D> coarseAdvect(coarseLevelSet, velocities);
      coarseAdvect.setAdvectionTime(advectionTime);
      
      // Use high-order schemes for accurate transport on coarse grid
      coarseAdvect.setSpatialScheme(SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
      coarseAdvect.setTemporalScheme(TemporalSchemeEnum::FORWARD_EULER);
      
      std::cout << "NeuralAdvect: Coarse points before Advect: " << coarseLevelSet->getNumberOfPoints() << std::endl;
      
      coarseAdvect.apply();

      {
        auto mesh = Mesh<T>::New();
        ToMesh<T, D>(coarseLevelSet, mesh).apply();
        VTKWriter<T>(mesh, "coarseLevelSet_afterAdvect.vtu").apply();
        ToSurfaceMesh<T, D>(coarseLevelSet, mesh).apply();
        VTKWriter<T>(mesh, "coarseLevelSet_afterAdvect.vtp").apply();
      }
      std::cout << "NeuralAdvect: Coarse points after Advect: " << coarseLevelSet->getNumberOfPoints() << std::endl;

      // Expand coarse level set to provide more context for the NN
      // and avoid "cliff" artifacts where the narrow band ends.
      Expand<T, D>(coarseLevelSet, 7).apply();

      // 3. NEURAL UPSAMPLE
      if (nnInference) {
          // If a NN is provided, use it to hallucinate details back to fine grid
          nnInference(coarseLevelSet, fineLevelSet);
      } else {
          // Fallback: Linear Interpolation via Resample
          // This will be blurry but valid
          Resample<T, D> upsampler(coarseLevelSet, fineLevelSet, fineDelta);
          upsampler.apply();
      }
      
      // Final cleanup
      fineLevelSet->finalize(2);
  }
};

} // namespace viennals
