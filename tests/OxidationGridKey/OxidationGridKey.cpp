// The per-node maps of the oxidation solvers (pressure, stress and
// concentration, the mask carry velocity) are keyed by
// detail::gridIndexHash alone, so it must be injective.  The hash_combine
// it replaced collided for 9 in 10 indices of a 5 nm trench domain, and
// nodes read each other's stress.
#include <lsOxidationSolverBase.hpp>
#include <lsTestAsserts.hpp>

#include <cstdio>
#include <unordered_set>

namespace ls = viennals;

// Every index of the cube [-extent, extent]^D gets its own key.
template <int D> void testInjective(int extent) {
  std::unordered_set<std::size_t> keys;
  std::size_t count = 0;
  viennahrle::Index<D> index;
  for (int d = 0; d < D; ++d)
    index[d] = -extent;
  while (true) {
    keys.insert(ls::detail::gridIndexHash<D>(index));
    ++count;
    int d = 0;
    for (; d < D; ++d) {
      if (index[d] < extent) {
        ++index[d];
        break;
      }
      index[d] = -extent;
    }
    if (d == D)
      break;
  }
  std::printf("%d-D: %zu indices, %zu keys\n", D, count, keys.size());
  VC_TEST_ASSERT(keys.size() == count)
}

viennahrle::Index<3> index3(viennahrle::IndexType x, viennahrle::IndexType y,
                            viennahrle::IndexType z) {
  viennahrle::Index<3> index;
  index[0] = x;
  index[1] = y;
  index[2] = z;
  return index;
}

int main() {
  testInjective<2>(400);
  testInjective<3>(60);
  // the ends of the packed range
  constexpr viennahrle::IndexType lo3 = -(1 << 20), hi3 = (1 << 20) - 1;
  VC_TEST_ASSERT(ls::detail::gridIndexHash<3>(index3(lo3, 0, 0)) !=
                 ls::detail::gridIndexHash<3>(index3(hi3, 0, 0)))
  VC_TEST_ASSERT(ls::detail::gridIndexHash<3>(index3(0, lo3, hi3)) !=
                 ls::detail::gridIndexHash<3>(index3(0, hi3, lo3)))
  return 0;
}
