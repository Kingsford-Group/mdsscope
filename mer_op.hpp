
#ifndef MER_OP_H
#define MER_OP_H

#include <inttypes.h>
#include <type_traits>
#include <cstdint>

typedef uint16_t mer_type;

template<typename T>
constexpr typename std::enable_if<std::is_integral<T>::value, T>::type
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

template<unsigned int k_, unsigned int alpha_>
struct mer_op_type {
    typedef mer_type mer_t;

    constexpr static unsigned int k = k_;
    constexpr static unsigned int alpha = alpha_;
    constexpr static mer_t nb_mers = ipow((mer_t)alpha, k);
    constexpr static mer_t nb_fmoves = nb_mers / alpha;

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

    static mer_t weight(const mer_t m) {
        mer_t w = 0;
        mer_t left = m;
        for(unsigned int i = 0; i < k; ++i, left /= alpha)
            w += left % alpha;
        return w;
    }
};

#endif // MER_OP_H
