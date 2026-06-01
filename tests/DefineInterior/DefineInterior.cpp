#include <lsDomain.hpp>
#include <lsInterior.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>
#include <hrleSparseStarIterator.hpp>
#include <iostream>

using namespace viennals;

template <int D> void testInteriorDefinition(const std::string &suffix,
                                              size_t expectedBefore,
                                              size_t expectedAfter) {
  using NumericType = double;

  auto circle = SmartPointer<Domain<NumericType, D>>::New(1.0);
  NumericType origin[D];
  for (int i = 0; i < D; ++i) origin[i] = 0.0;
  auto sphere = SmartPointer<Sphere<NumericType, D>>::New(origin, 5.0);
  MakeGeometry<NumericType, D>(circle, sphere).apply();

  auto &domain = circle->getDomain();
  auto &grid = circle->getGrid();
  viennahrle::ConstSparseStarIterator<typename Domain<NumericType, D>::DomainType, 1>
      itBefore(domain, grid.getMinGridPoint());
  viennahrle::Index<D> endIdx = grid.incrementIndices(grid.getMaxGridPoint());
  size_t interiorBefore = 0;
  while (itBefore.getIndices() < endIdx) {
    if (itBefore.getCenter().isDefined() && itBefore.getCenter().getValue() < 0) {
      interiorBefore++;
    }
    itBefore.next();
  }

  VC_TEST_ASSERT(interiorBefore == expectedBefore)

  {
    auto mesh = SmartPointer<Mesh<NumericType>>::New();
    ToMesh<NumericType, D>(circle, mesh).apply();
    VTKWriter<NumericType>(mesh, "interior_before_" + suffix + ".vtk").apply();
  }

  Interior<NumericType, D> filling(circle);
  filling.apply();

  LSTEST_ASSERT_VALID_LS(circle, NumericType, D)

  viennahrle::ConstSparseStarIterator<typename Domain<NumericType, D>::DomainType, 1>
      itAfter(circle->getDomain(), grid.getMinGridPoint());
  size_t interiorAfter = 0;
  while (itAfter.getIndices() < endIdx) {
    if (itAfter.getCenter().isDefined() && itAfter.getCenter().getValue() < 0) {
      interiorAfter++;
    }
    itAfter.next();
  }

  {
    auto mesh = SmartPointer<Mesh<NumericType>>::New();
    ToMesh<NumericType, D>(circle, mesh).apply();
    VTKWriter<NumericType>(mesh, "interior_after_" + suffix + ".vtk").apply();
  }

  VC_TEST_ASSERT(interiorAfter == expectedAfter)

  std::cout << "D=" << D << ": " << interiorBefore << " → " << interiorAfter
            << " interior points\n";
}

int main() {
  testInteriorDefinition<2>("2D", 24, 69);
  testInteriorDefinition<3>("3D", 210, 485);

  return 0;
}
