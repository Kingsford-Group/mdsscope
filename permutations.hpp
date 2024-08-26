#ifndef PERMUTATIONS_H_
#define PERMUTATIONS_H_

#include "mer_op.hpp"
#include <limits>
#include <random>
#include <cstdint>
#include <type_traits>
#include <xxhash.h>

// xxhash function with a seed. Computed as a uint64_t word, then folded down
// (XOR high part of word onto lower parts) to type T if needed.
template<typename T>
struct xxhash {
  static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value, "Not an unsigned integer type");
  const uint64_t seed;

  template<typename RNG>
  xxhash(RNG& rng)
  : seed(std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(rng))
  {}

  T operator()(const T& x) const {
    uint64_t h = XXH64((const void*)&x, sizeof(x), seed);

    // Fold the word on itself until sizeof(T)
    for(unsigned i = sizeof(uint64_t) * 8; i > sizeof(T) * 8; i /= 2) {
      h ^= h >> i;
    }

    return h;
  }
};

template<typename T>
struct FeistelPermutation {
  const xxhash<T> hash;

  template<typename RNG>
  FeistelPermutation(RNG& rng)
  : hash(rng)
  {}

  T operator()(const T& x) const {
    constexpr unsigned wbit = sizeof(T) * 8;
    constexpr unsigned hbit = wbit / 2;
    constexpr T lhalf = (T)-1 >> hbit;
    // constexpr T hhalf = lhalf << hbit;

    const T h = hash(x & lhalf);
    return ((x & lhalf) << hbit) | ((h & lhalf) ^ (h >> hbit) ^ (x >> hbit));
  }
};

template<typename T>
struct LubyRackofPermutation {
  const FeistelPermutation<T> DF1, DF2, DF3, DF4;

  template<typename RNG>
  LubyRackofPermutation(RNG& rng)
  : DF1(rng)
  , DF2(rng)
  , DF3(rng)
  , DF4(rng)
  {}

  T operator()(const T& x) const {
    return DF4(DF3(DF2(DF1(x))));
  }
};

#endif // PERMUTATIONS_H_
