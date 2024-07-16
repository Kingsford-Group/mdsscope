
#ifndef MER_OP_H
#define MER_OP_H

// #include <inttypes.h>
#include <type_traits>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <numeric>
#include <ostream>
#include <iostream>

// Print 128 bit long integers. Not very fast. Ignore formatting
std::ostream& operator<<(std::ostream& os, __uint128_t x) {
	static constexpr int buflen = 40;
	char buf[buflen];
	char* ptr = &buf[buflen - 1];
	*ptr = '\0';

    do {
		--ptr;
        *ptr = ((char)(x % 10 + '0'));
        x /= 10;
    } while(x > 0);
	return os << ptr;
}

// Number of bits to encode a^k
constexpr unsigned int log2ak(unsigned int a, unsigned k) {
    std::size_t val = 1;
    unsigned int res = 0;
    std::size_t bound = std::numeric_limits<std::size_t>::max() / a;
    while(k > 0) {
        if(val < bound) {
            val *= a;
            k -= 1;
        } else {
            val /= 2;
            res += 1;
        }
    }
    while(val > 0) {
        val /= 2;
        res += 1;
    }
    return res;
}

// Find the integer type in the list T, Ts... that is large enough for this number of bits
template<unsigned bits, typename T, typename... Ts>
struct optimal_int {
    static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value, "Not an unsigned integer type");
    typedef typename std::conditional<(bits > 8 * sizeof(T)), typename optimal_int<bits, Ts...>::type, T>::type type;
};

template<unsigned bits, typename T>
struct optimal_int<bits, T> {
    static_assert(bits <= 8 * sizeof(T), "Too many bits");
    typedef T type;
};

template<typename T>
// constexpr typename std::enable_if<std::is_integral<T>::value, T>::type
constexpr T
ipow(T base, unsigned int exp) {
    T result = 1;
    for(;;) {
        if(exp & 1) result *= base;
        exp >>= 1;
        if(!exp) break;
        base *= base;
    }

    return result;
}

constexpr size_t nb_necklaces(unsigned a, unsigned k) {
	size_t res = 0;
    for(unsigned i = 1; i <= k; ++i) {
		res += ipow(a, std::gcd(i, k));
	}
	return res / k;
}

// Bit twiddling operation to reverse bits in a word. Used for reverse
// complementation when the alphabet is a power of 2 (alpha == 2 or 4). Only
// defined on unsigned word types.
//
// Checkered mask. cmask<uint16_t, 1> is every other bit on (0x55).
// cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
struct cmask {
  static constexpr
  std::enable_if<std::is_unsigned<U>::value && std::is_integral<U>::value, U>::type
  v = (cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
};

// When len is half of the word size, shifting by (2 * len) is undefined
// behavior. Fix it here.
template<typename U, int l>
struct cmask<U, sizeof(U) * 4, l> {
  static constexpr
  std::enable_if<std::is_unsigned<U>::value && std::is_integral<U>::value, U>::type
  v = (((U)1 << (sizeof(U) * 4)) - 1);
};

// Base case, when l = 0, start with empty (0) word.
template<typename U, int len>
struct cmask<U, len, 0> {
  static constexpr
  std::enable_if<std::is_unsigned<U>::value && std::is_integral<U>::value, U>::type
  v = 0;
};

template<typename U, unsigned alpha>
inline
std::enable_if<std::is_unsigned<U>::value && std::is_integral<U>::value, U>::type
word_reverse_complement(U w) {
    if constexpr (alpha == 2)
        w = ((w >> 1)  & cmask<U, 1 >::v) | ((w & cmask<U, 1 >::v) << 1);
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    if constexpr (sizeof(U) >= 2)
        w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    if constexpr (sizeof(U) >= 4)
        w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    if constexpr (sizeof(U) >= 8)
        w = ((w >> 32) & cmask<U, 32>::v) | ((w & cmask<U, 32>::v) << 32);
    if constexpr (sizeof(U) >= 16)
        w = ((w >> 64) & cmask<U, 64>::v) | ((w & cmask<U, 64>::v) << 64);
    // return ~w;
  return ((U)-1) - w;
}

template<unsigned int k_, unsigned int alpha_>
struct mer_op_type {
//    typedef mer_type mer_t;
    constexpr static unsigned int ak_bits = log2ak(alpha_, k_);
    typedef typename optimal_int<ak_bits, uint8_t, uint16_t, uint32_t, uint64_t, __uint128_t>::type mer_t;

    // Skip many program if encoding k takes too many bits
    constexpr static unsigned int max_bits = 34;

    constexpr static unsigned int k = k_;
    constexpr static unsigned int alpha = alpha_;
    constexpr static mer_t nb_mers = ipow((mer_t)alpha, k);
    constexpr static mer_t nb_fmoves = nb_mers / alpha;
    constexpr static mer_t nb_necklaces = ::nb_necklaces(alpha, k);

    // Left base
    static inline mer_t lb(const mer_t m) {
        return m / nb_fmoves;
    }

    // Right base
    static inline mer_t rb(const mer_t m) {
        return m % alpha;
    }

    // Next mer by pure rotation
    static inline mer_t nmer(const mer_t m) {
        return nmer(m, lb(m));
    }

    // Next mer using given new base (do an F-move if base != lb(m))
    static inline mer_t nmer(const mer_t m, const mer_t base) {
        return ((m * alpha) % nb_mers) + (base % alpha);
    }

    // Previous mer by pure rotation
    static inline mer_t pmer(const mer_t m) {
        return pmer(m, rb(m));
    }

    // Previous mer given new base (do a RF-move if base != lr(m))
    static inline mer_t pmer(const mer_t m, const mer_t base) {
        return (m / alpha) + (base % alpha) * nb_fmoves;
    }

    // Fmove corresponding to mer
    static inline mer_t fmove(const mer_t m) {
        return m % nb_fmoves;
    }

    // RFmove corresponding to mer
    static inline mer_t rfmove(const mer_t m) {
        return m / alpha;
    }

    // Are two mers left companions?
    static inline bool are_lc(const mer_t m1, const mer_t m2) {
        return fmove(m1) == fmove(m2);
    }

    // Are two mers right companions?
    static inline bool are_rc(const mer_t m1, const mer_t m2) {
        return (m1 / alpha) == (m2 / alpha);
    }

    // Left-companion with that base
    static inline mer_t lc(const mer_t m, const mer_t base) {
        return fmove(m) + (base % alpha) * nb_fmoves;
    }

    // Right-companion with that base
    static inline mer_t rc(const mer_t m, const mer_t base) {
        return m - rb(m) + (base % alpha);
    }

    static mer_t homopolymer(const mer_t base) {
        mer_t m = 0;
        const auto b = base % alpha;
        for(unsigned int i = 0; i < k; ++i)
            m = m * alpha + b;
        return m;
    }

    static mer_t is_homopolymer(const mer_t m) {
        // Same as period-1 PCR
        return nmer(m) == m;
    }

    static mer_t is_homopolymer_fm(const mer_t fm) {
        return ((fm * alpha) % nb_fmoves + (fm % alpha)) == fm;
    }

    static mer_t weight(const mer_t m) {
        mer_t w = 0;
        mer_t left = m;
        for(unsigned int i = 0; i < k; ++i, left /= alpha)
            w += left % alpha;
        return w;
    }

    inline static mer_t reverse_comp(const mer_t m) {
        // Optimize for binary and DNA alphabet with bit twiddling
        if constexpr (alpha == 2 || alpha == 4) {
            const mer_t wc = word_reverse_complement<mer_t,alpha>(m);
            constexpr unsigned shift = (8 * sizeof(mer_t) - (alpha/2)*k);
            return  wc >> shift;
        } else {
            // For general alphabets, looping algorithm
            mer_t res = 0;
            mer_t left = m;
            for(unsigned int i = 0; i < k; ++i, left /= alpha)
                res = (res * alpha) + (alpha - 1 - (left % alpha));
            return res;
        }
    }

    static mer_t canonical(const mer_t m) {
        return std::min(m, reverse_comp(m));
    }
};


template<unsigned int k_, unsigned int alpha_>
struct amer_type {
    typedef mer_op_type<k_, alpha_> mer_ops;
    typedef mer_ops::mer_t mer_t;
    mer_t val;

    amer_type() = default;
    amer_type(mer_t x) : val(x) {}
    amer_type(const amer_type& rhs) : val(rhs.val) {}
    static amer_type homopolymer(const mer_t base) { return mer_ops::homopolymer(base); }

    inline amer_type lb() const { return mer_ops::lb(val); }
    inline amer_type rb() const { return mer_ops::rb(val); }
    inline amer_type nmer() const { return mer_ops::nmer(val); }
    inline amer_type nmer(const mer_t base) const { return mer_ops::nmer(val, base); }
    inline amer_type pmer() const { return mer_ops::pmer(val); }
    inline amer_type pmer(const mer_t base) const { return mer_ops::pmer(val, base); }
    inline amer_type fmove() const { return mer_ops::fmove(val); }
    inline amer_type rfmove() const { return mer_ops::rfmove(val); }
    inline bool are_lc(const amer_type& rhs) const { return mer_ops::are_lc(val, rhs.val); }
    inline bool ar_rc(const amer_type& rhs) const { return mer_ops::are_rc(val, rhs.val); }
    inline amer_type lc(const mer_t base) const { return mer_ops::lc(val, base); }
    inline amer_type rc(const mer_t base) const { return mer_ops::rc(val, base); }
    inline bool is_homopolymer() const { return mer_ops::is_homopolymer(val); }
    inline bool is_homopolymer_fm() const { return mer_ops::is_homopolymer_fm(val); }
    inline mer_t weight() const { return mer_ops::weight(val); }
    inline amer_type reverse_comp() const { return mer_ops::reverse_comp(val); }
    inline amer_type canonical() const { return mer_ops::canonical(val); }

    struct mer_rc_pair {
        amer_type mer, rc;
        mer_rc_pair(const amer_type& m)
            : mer(m)
            , rc(m.reverse_comp())
            {}
        mer_rc_pair& operator++() {
            ++mer.val;
            rc.val -= mer_ops::nb_fmoves;
            return *this;
        }
        mer_rc_pair& operator--() {
            --mer;
            rc += mer_ops::nb_fmoves;
            return *this;
        }
    };

};

template<unsigned int k_, unsigned int alpha_>
std::ostream& operator<<(std::ostream& os, const amer_type<k_, alpha_>& m) {
    return os << (uint64_t)m.val;
}

template<unsigned int k_, unsigned int alpha_>
inline bool operator==(const amer_type<k_, alpha_>& x, const amer_type<k_, alpha_>& y) {
    return x.val == y.val;
}

template<unsigned int k_, unsigned int alpha_>
inline bool operator!=(const amer_type<k_, alpha_>& x, const amer_type<k_, alpha_>& y) {
    return x.val != y.val;
}

template<unsigned int k_, unsigned int alpha_>
inline bool operator<(const amer_type<k_, alpha_>& x, const amer_type<k_, alpha_>& y) {
    return x.val < y.val;
}


template<unsigned int k_, unsigned int alpha_>
struct std::hash<amer_type<k_, alpha_>> {
    typedef amer_type<k_, alpha_> amer_t;
    std::size_t operator()(const amer_t& x) const noexcept {
        return std::hash<typename amer_t::mer_t>{}(x.val);
    }
};

template<unsigned int k_, unsigned int alpha_>
inline std::ostream& operator<<(std::ostream& os, const typename amer_type<k_, alpha_>::mer_rc_pair& pair) {
    return os << '(' << pair.mer << ',' << pair.rc << ')';
}

#endif // MER_OP_H
