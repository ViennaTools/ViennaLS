#ifndef LS_CONCEPTS_HPP
#define LS_CONCEPTS_HPP

#include <cstddef>
#include <type_traits>

namespace lsConcepts {
// use any type that can be assigned any value (so anything but void)
using AssignType = std::nullptr_t;
// some value that can be used as the default parameter
inline constexpr AssignType assignable = AssignType();

template <class Base, class Derived>
using IsBaseOf =
    std::enable_if_t<std::is_base_of<Base, Derived>::value, AssignType>;

template <class A, class B>
using IsSame = std::enable_if_t<std::is_same<A, B>::value, AssignType>;

template <class A, class B>
using IsNotSame = std::enable_if_t<!std::is_same<A, B>::value, AssignType>;

template <class T>
using IsFloatingPoint =
    std::enable_if_t<std::is_floating_point<T>::value, AssignType>;

} // namespace lsConcepts

#endif // LS_CONCEPTS_HPP