#ifndef MER_OP_H
#define MER_OP_H

#include <inttypes.h>
#include <type_traits>

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
    constexpr static mer_type nb_mers = ipow((mer_type)alpha, k);
    constexpr static mer_type nb_fmoves = nb_mers / alpha;

    // Left base
    static inline mer_type lb(const mer_type m) {
        return m / nb_fmoves;
    }

    // Right base
    static inline mer_type rb(const mer_type m) {
        return m % alpha;
    }

    // Next mer by pure rotation
    static inline mer_type nmer(const mer_type m) {
        return nmer(m, lb(m));
    }

    // Next mer using given new base (do an F-move if base != lb(m))
    static inline mer_type nmer(const mer_type m, const mer_type base) {
        return ((m * alpha) % nb_mers) + (base % alpha);
    }

    // Previous mer by pure rotation
    static inline mer_type pmer(const mer_type m) {
        return pmer(m, rb(m));
    }

    // Previous mer given new base (do a RF-move if base != lr(m))
    static inline mer_type pmer(const mer_type m, const mer_type base) {
        return (m / alpha) + (base % alpha) * nb_fmoves;
    }

    // Fmove corresponding to mer
    static inline mer_type fmove(const mer_type m) {
        return m % nb_fmoves;
    }

    // RFmove corresponding to mer
    static inline mer_type rfmove(const mer_type m) {
        return m / alpha;
    }

    // Are two mers left companions?
    static inline bool are_lc(const mer_type m1, const mer_type m2) {
        return fmove(m1) == fmove(m2);
    }

    // Are two mers right companions?
    static inline bool are_rc(const mer_type m1, const mer_type m2) {
        return (m1 / alpha) == (m2 / alpha);
    }

    // Left-companion with that base
    static inline mer_type lc(const mer_type m, const mer_type base) {
        return fmove(m) + (base % alpha) * nb_fmoves;
    }

    // Right-companion with that base
    static inline mer_type rc(const mer_type m, const mer_type base) {
        return m - rb(m) + (base % alpha);
    }

    static mer_type homopolymer(const mer_type base) {
        mer_type m = 0;
        const auto b = base % alpha;
        for(unsigned int i = 0; i < k; ++i)
            m = m * alpha + b;
        return m;
    }

    static mer_type weight(const mer_type m) {
        mer_type w = 0;
        mer_type left = m;
        for(unsigned int i = 0; i < k; ++i, left /= alpha)
            w += left % alpha;
        return w;
    }
};

#endif // MER_OP_H
