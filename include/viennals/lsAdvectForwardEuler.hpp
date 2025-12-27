#pragma once

#include <lsAdvect.hpp>

namespace viennals {

template <class T, int D> class AdvectForwardEuler : public Advect<T, D> {
public:
  using Advect<T, D>::Advect;
};

PRECOMPILE_PRECISION_DIMENSION(AdvectForwardEuler);

} // namespace viennals
