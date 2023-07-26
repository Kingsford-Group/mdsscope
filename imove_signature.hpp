#ifndef IMOVE_SIGNATURE_H_
#define IMOVE_SIGNATURE_H_

#include <vector>
#include <unordered_map>
#include <xxh3.h>
#include <ostream>
#include <istream>
#include <cstring>

#include "mer_op.hpp"
#include "common.hpp"

// Base class that zeros itself upon construction. Guarantees that even the
// unused bits (e.g., because of alignment of the fields) are deterministically
// set to 0. Makes the hash function below work.
template<typename D>
struct zeroed {
    zeroed() {
        memset(this, '\0', sizeof(D));
    }
};

template<typename mer_op_type>
struct imove_type : zeroed<imove_type<mer_op_type>> {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;
    typedef uint8_t mask_type;

    static constexpr mask_type set_mask(mask_type x, unsigned int b) {
        return b ? set_mask(x | ((mask_type)1 << (b - 1)), b-1) : x;
    }

    // TODO: interpretation of mask?
    static constexpr mask_type none = set_mask(0, mer_op_t::alpha);
    static constexpr mask_type all = none;

    // An I-move is identified by the F-move suffix/prefix. im a a bit mask of
    // length ALPHA. A bit of 0 means at position x the left-companion lc(fm, x)
    // is in the set, while a bit of 1 means rc(fm, x) is in the set. A mask of
    // 0 does not represent an I-move but the F-move fm itself. A mask of all
    // ones (up to ALPHA) represents the reverse F-move fm. Any other mask is
    // one of the 2^ALPHA-2 valid I-move.
    //
    // Not that if an I-move is possible, then any submask is also a possible
    // I-move.
    mer_t fm; // F-move part
    mask_type im; // I-move part

    imove_type() : fm(-1), im(-1) { } // Invalid will need to update
    imove_type(mer_t _fm, mask_type _im) : fm(_fm), im(_im) {}
};

template<typename mer_op_type>
bool operator==(const imove_type<mer_op_type>& im1, const imove_type<mer_op_type>& im2) {
    return im1.fm == im2.fm && im1.im == im2.im;
}

template<typename mer_op_type>
using imove_sig_type = std::vector<imove_type<mer_op_type>>;

template<typename mer_op_type>
std::ostream& operator<<(std::ostream& os, const imove_type<mer_op_type>& rhs) {
    return os << rhs.fm << ':' << (unsigned int)rhs.im;
}

template<typename mer_op_type>
std::istream& operator>>(std::istream& is, imove_type<mer_op_type>& rhs) {
    typename mer_op_type::mer_t fm;
    char sep;
    unsigned int im;
    is >> fm >> sep >> im;
    if(fm < mer_op_type::nb_fmoves && im < imove_type<mer_op_type>::none && sep == ':') {
        rhs.fm = fm;
        rhs.im = im;
    } else {
        is.setstate(std::ios::failbit);
    }
    return is;
}

namespace std {
template<typename mer_op_type>
struct hash<imove_sig_type<mer_op_type>> {
    static constexpr uint64_t seed = 0xd336ea32c33ce21fUL;
    size_t operator()(const imove_sig_type<mer_op_type>& sig) const {
        return XXH64(sig.data(), sig.size() * sizeof(imove_type<mer_op_type>), seed);
    }
};
}

template<typename mer_op_type>
using signatures_type = std::unordered_map<imove_sig_type<mer_op_type>, size_t>;

#endif // IMOVE_SIGNATURE_H_
