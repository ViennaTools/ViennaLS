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

  typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  return ls::SmartPointer<ls::Domain<double, 3>>::New(bounds, boundaryCons,
                                                      gridDelta);
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
