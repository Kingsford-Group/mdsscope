#ifndef DBG_H_
#define DBG_H_

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define E(e) (TOSTRING(e) "<") << (e) << '>'

// Own assert with message
#ifdef NDEBUG
#define assert2(condition, message) do { } while(false)
#else
#include <exception>
#include <stdexcept>
#include <iostream>

#define assert2(condition, message)                                     \
    do {                                                                \
    if(!(condition)) {                                                  \
    std::cerr << "Assertion failed: (" << #condition << "), "           \
    << __FUNCTION__                                                     \
    << '@' << __FILE__                                                  \
    << ':' << __LINE__ << '.'                                           \
    << '\n' << message << std::endl;                                    \
    throw std::logic_error("Assert2");                                  \
    }                                                                   \
     } while(false)
#endif

#endif // DBG_H_
