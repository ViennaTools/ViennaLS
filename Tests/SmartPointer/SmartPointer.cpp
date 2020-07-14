#include <iostream>
#include <lsDomain.hpp>
#include <lsSmartPointer.hpp>

constexpr int D = 3;
typedef double NumericType;

auto makeLSDomain() {
  double gridDelta = 1.1;

  double extent = 50;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  return lsSmartPointer<lsDomain<double, 3>>::New(bounds, boundaryCons,
                                                  gridDelta);
}

unsigned long
getNumberOfPoints(lsSmartPointer<const lsDomain<NumericType, D>> domain) {
  return domain->getNumberOfPoints();
}

int main() {
  auto domain = makeLSDomain();

  // lsSmartPointer<const lsDomain<NumericType, D>> constDomain = domain;

  std::cout << getNumberOfPoints(domain) << std::endl;

  std::cout << "Number of Points: " << domain->getNumberOfPoints() << std::endl;
  std::cout << "Number of references to lsDomain: " << domain.use_count()
            << std::endl;
  {
    auto domain2 = domain;
    std::cout << "Number of references to lsDomain: " << domain.use_count()
              << std::endl;
  }
  std::cout << "Number of references to lsDomain: " << domain.use_count()
            << std::endl;

  return 0;
}
