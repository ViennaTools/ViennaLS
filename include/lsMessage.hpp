#ifndef LS_MESSAGE_HPP
#define LS_MESSAGE_HPP

#include <iostream>

/// Singleton class for thread-safe logging.
class lsMessage {
  std::string message;

  bool error = false;
  const unsigned tabWidth = 4;

  lsMessage() {}

public:
  // delete constructors to result in better error messages by compilers
  lsMessage(const lsMessage &) = delete;
  void operator=(const lsMessage &) = delete;

  static lsMessage &getInstance() {
    static lsMessage instance;
    return instance;
  }

  lsMessage &add(std::string s) {
#pragma omp critical
    { message += "\n" + std::string(tabWidth, ' ') + "WARNING: " + s + "\n"; }
    return *this;
  }

  lsMessage &addWarning(std::string s) {
#pragma omp critical
    { message += "\n" + std::string(tabWidth, ' ') + "WARNING: " + s + "\n"; }
    return *this;
  }

  lsMessage &addError(std::string s, bool shouldAbort = true) {
#pragma omp critical
    { message += "\n" + std::string(tabWidth, ' ') + "ERROR: " + s + "\n"; }
    // always abort once error message should be printed
    error = true;
    // abort now if asked
    if (shouldAbort)
      print();
    return *this;
  }

  lsMessage &addDebug(std::string s) {
#pragma omp critical
    { message += std::string(tabWidth, ' ') + "DEBUG: " + s + "\n"; }
    return *this;
  }

  void print(std::ostream &out = std::cout) {
    out << message;
    message.clear();
    if (error)
      abort();
  }
};

#endif // LS_MESSAGE_HPP
