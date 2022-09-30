#ifndef LS_SMART_POINTER_HPP
#define LS_SMART_POINTER_HPP

#include <memory>

/// std::shared_ptr wrapper for use with ViennaLS.
/// lsSmartPointers should be created using the function ::New(...).
/// All other interface functions are identical to std::shared_ptr
template <class T> class lsSmartPointer : public std::shared_ptr<T> {
public:
  // lsSmartPointer(T& passedObject) :
  // std::shared_ptr<T>(std::make_shared(passedObject)) {}

  // Make visible all constructors of std::shared_ptr
  // including copy constructors
  template <typename... Args>
  lsSmartPointer(Args &&...args)
      : std::shared_ptr<T>(std::forward<Args>(args)...) {}

  /// Use this function to create new objects when using ViennaLS
  template <typename... TArgs> static lsSmartPointer New(TArgs &&...targs) {
    return lsSmartPointer(std::make_shared<T>(std::forward<TArgs>(targs)...));
  }
};

#endif // LS_SMART_POINTER_HPP