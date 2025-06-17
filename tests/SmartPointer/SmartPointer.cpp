#include <iostream>
#include <lsDomain.hpp>

namespace ls = viennals;

constexpr int D = 3;
typedef double NumericType;

auto makeLSDomain() {
  double gridDelta = 1.1;

  double extent = 50;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename ls::BoundaryConditionEnum boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

  return ls::Domain<double, 3>::New(bounds, boundaryCons, gridDelta);
}

unsigned long
getNumberOfPoints(ls::SmartPointer<const ls::Domain<NumericType, D>> domain) {
  return domain->getNumberOfPoints();
}

int main() {
  auto domain = makeLSDomain();

  // ls::SmartPointer<const ls::Domain<NumericType, D>> constDomain = domain;

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
