#ifndef LS_CONCEPTS_HPP
#define LS_CONCEPTS_HPP

#include <cstddef>

namespace lsConcepts {
// use any type that can be assigned any value (so anything but void)
using AssignType = std::nullptr_t;
// some value that can be assigned to nullptr_t
constexpr AssignType assignable{};

template <class Base, class Derived>
using IsBaseOf = typename std::enable_if<std::is_base_of<Base, Derived>::value,
                                         AssignType>::type;

template <class A, class B>
using IsSame =
    typename std::enable_if<std::is_same<A, B>::value, AssignType>::type;

template <class A, class B>
using IsNotSame =
    typename std::enable_if<!std::is_same<A, B>::value, AssignType>::type;

} // namespace lsConcepts

#endif // LS_CONCEPTS_HPP