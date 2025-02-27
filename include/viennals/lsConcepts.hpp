#pragma once

#include <cstddef>
#include <type_traits>

namespace lsConcepts {
// use any type that can be assigned any value (so anything but void)
using AssignType = std::nullptr_t;
// some value that can be used as the default parameter
inline constexpr AssignType assignable = AssignType();

template <class Base, class Derived>
using IsBaseOf = std::enable_if_t<std::is_base_of_v<Base, Derived>, AssignType>;

template <class A, class B>
using IsSame = std::enable_if_t<std::is_same_v<A, B>, AssignType>;

template <class A, class B>
using IsNotSame = std::enable_if_t<!std::is_same_v<A, B>, AssignType>;

template <class T>
using IsFloatingPoint =
    std::enable_if_t<std::is_floating_point_v<T>, AssignType>;

} // namespace lsConcepts
