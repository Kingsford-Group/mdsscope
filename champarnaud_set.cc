#include <cstdlib>
#include <bitset>
#include <tuple>
#include <limits>
#include <algorithm>
#include <vector>

#include "champarnaud_set.hpp"
#include "common.hpp"
#include "dbg.hpp"


#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

// lyndon_test: Test if a work/k-mer is a Lyndon word. A word is a Lyndon word
// if it is primitive (not equal to a proper power: u^n with n > 1) and minimal
// (least in its conjugacy class / PCR).
//
// Equivalently, a word is a Lyndon word if it is strictly smaller than all its
// proper non-empty suffixes. This is the test used here. Use template recursive
// expression for the compiler to optimize the modulo (%) operation.
template<unsigned i, mer_t pmask, mer_t smask>
struct lyndon_test {
    constexpr static bool is_lyndon(const mer_t m) {
        // std::cout << "is_lyndon " << i << ' ' << (size_t)mask << ' ' << (size_t)m << ' ' << (size_t)(m % mask) << ' ' << (size_t)suffix << '\n';
        // Word read MSB to LSB. Division by alpha^{i} gives the i-th prefix, modulo alpha^{k-i} gives the i-th suffix
        return ((m / pmask) < (m % smask)) && lyndon_test<i - 1, pmask * mer_ops::alpha, smask / mer_ops::alpha>::is_lyndon(m);
    }
};

template<mer_t pmask, mer_t smask>
struct lyndon_test<0, pmask, smask> {
    constexpr static bool is_lyndon(const mer_t m) {
        // static_assert(pmask == mer_ops::nb_mers);
        // std::cout << "is_lyndon " << 1 << ' ' << (size_t)mask << ' ' << (size_t)m << ' ' << (size_t)suffix << '\n';
        return true;
    }
};

template<unsigned wlen = mer_ops::k>
inline constexpr bool is_lyndon(mer_t m) { return lyndon_test<wlen - 1, mer_ops::alpha, ipow(mer_ops::alpha, wlen - 1)>::is_lyndon(m); }

// Find largest n such that w = l^n u, where l is a Lyndon word and for some
// word u. oi is the length of the Lyndon word and i = n*oi. l is the lyndon
// word and ln is l^n.
template<unsigned i, unsigned oi, mer_t mask, mer_t omask, bool underk>
struct word_divide { };

template<unsigned i, unsigned oi, mer_t mask, mer_t mstep>
struct word_divide<i, oi, mask, mstep, true> {
    constexpr static std::tuple<mer_t, unsigned> divide(mer_t w, mer_t l, mer_t ln) {
        constexpr mer_t nmask = mask / mstep;
        const mer_t lnp1 = ln * mstep + l;
        // std::cout << "\tword_divide under " << i << ' ' << oi << ' ' << (size_t)mask << ' ' << (size_t)mstep << ' ' << (size_t)l << ' ' << (size_t)ln << '\n'
        //     << "\t\t" << (size_t)nmask << ' ' << (size_t)w << ' ' << (size_t)(w / nmask) << ' ' << (size_t)lnp1 << '\n';
        if(w / nmask == lnp1)
            return word_divide<i + oi, oi, nmask, mstep, (i + oi <= mer_ops::k)>::divide(w, l, lnp1);
        constexpr mer_t smask = mer_ops::nb_mers / mask; // suffix mask
        return std::make_tuple(ln + (w % mask) * smask, i - oi);
    }
};

template<unsigned i, unsigned oi, mer_t mask, mer_t mstep>
struct word_divide<i, oi, mask, mstep, false> {
    constexpr static std::tuple<mer_t, unsigned> divide(mer_t w, mer_t l, mer_t ln) {
        // std::cout << "\tword_divide over " << i << ' ' << oi << ' ' << (size_t)mask << ' ' << (size_t)mstep << ' ' << (size_t)l << ' ' << (size_t)ln << '\n';
        constexpr mer_t smask = mer_ops::nb_mers / mask;
        return std::make_tuple(ln + (w % mask) * smask, i - oi);
    }
};

template<unsigned i, mer_t mask>
struct primary_division {
    constexpr static mer_t divide(mer_t w) {
        const mer_t l = w / mask;
        if(is_lyndon<i>(l)) {
            // Prefix of length i is a Lyndon word. Do division and return value
            // if length of rest is less than length of Lyndon word prefix
            const auto res = word_divide<i + i, i, mask, mer_ops::nb_mers / mask, (i + i <= mer_ops::k)>::divide(w, l, l);
            assert2(mer_ops::k >= std::get<1>(res), "Length of primary part longet than k");
            const unsigned len_remain = mer_ops::k - std::get<1>(res);
            // std::cout << "\tlyndon prefix " << i << ' ' << (size_t)mask << ' ' << (size_t)w << ' ' << (size_t)l << ' ' << len_remain << '\n';
            if(len_remain < i) return std::get<0>(res);
        }
        return primary_division<i + 1, mask / mer_ops::alpha>::divide(w);
    }
};

template<mer_t mask>
struct primary_division<mer_ops::k + 1, mask> {
    // Error. Should not get there. Stop recursion and generate error.
    constexpr static mer_t divide(mer_t w) {
        assert2(false, "Primary_division failed, should not get to length k+1");
        return mer_ops::nb_mers;
    }
};

// Given a minimal mer (smallest mer in a PCR), return the champarnaud mer for
// that PCR.
mer_t champarnaud_mer(mer_t m) {
    return primary_division<1, mer_ops::nb_mers / mer_ops::alpha>::divide(m);
}

int main(int argc, char* argv[]) {
    champarnaud_set args(argc, argv);

    // std::cout << "4 5 " << is_lyndon<4>(5) << '\n';

    std::bitset<mer_ops::nb_mers> done;
    std::vector<mer_t> res;

    for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
        if(done.test(m)) continue; // Already done that PCR

        // By definition, m is the smallest. Print the PCR and whether m is a
        // Lyndon word.
        // std::cout << (size_t)m;
        done.set(m);
        for(mer_t nm = mer_ops::nmer(m); nm != m; nm = mer_ops::nmer(nm)) {
            // std::cout << ',' << (size_t)nm;
            done.set(nm);
        }
        // std::cout << ": " << is_lyndon(m) << '\n';
        const auto mds_member = champarnaud_mer(m);
        // std::cout << "->" << (size_t)mds_member << "<-\n";
        res.push_back(mds_member);
    }
    std::sort(res.begin(), res.end());
    std::cout << joinT<size_t>(res, ',') << '\n';

    return EXIT_SUCCESS;
}
