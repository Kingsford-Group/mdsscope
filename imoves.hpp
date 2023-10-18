#ifndef _IMOVES_H_
#define _IMOVES_H_

#include <vector>
#include <algorithm>
#include <bitset>
#include <cstdint>

#include "common.hpp"
#include "mer_op.hpp"
#include "mds_op.hpp"
#include "imove_signature.hpp"

// Encodes the constrains between lc/rc around an F-move.
template<typename mer_op_type>
struct f_constrains {
    typedef typename mer_op_type::mer_t mer_t;
    static constexpr size_t len = mer_op_type::alpha * (mer_op_type::alpha - 1);
    std::bitset<len> bits;

    inline static size_t index(mer_t i, mer_t j) { return i * (mer_op_type::alpha - 1) + j - (i < j); }
    bool get(mer_t i, mer_t j) const { return bits[index(i, j)]; }
    void set(mer_t i, mer_t j) { bits.set(index(i, j)); }
    void reset() { bits.reset(); }
};

template<typename mer_op_type>
std::ostream& operator<<(std::ostream& os, const f_constrains<mer_op_type>& c) {
    for(size_t i = 0; i < mer_op_type::alpha; ++i) {
        for(size_t j = 0; j < mer_op_type::alpha; ++j)
            os << (i != j ? (c.get(i, j) ? '1' : '0') : '.');
        os << '|';
    }
    return os;
}

// Figure out the possible I-moves by doing DFSs. Caches data: not multi-thread
// safe (should have 1 object per-thread).
template<typename mer_op_type>
struct imoves_type {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_t::mer_t mer_t;
    typedef imove_type<mer_op_type> imove_t;
    typedef imove_sig_type<mer_op_type> imove_sig_t;
    typedef typename imove_t::mask_type mask_t;

    std::vector<f_constrains<mer_op_type>> constrained;
    std::vector<tristate_t> visited;
    // std::vector<mer_t> mds;
    std::vector<bool> done;

    imoves_type()
        : constrained(mer_op_t::nb_fmoves)
        , visited(mer_op_t::nb_mers)
        // , mds(mer_op_t::nb_mers)
        , done(mer_op_t::nb_fmoves)
    { }

    // Fill the constrained vector marking the edges (potential i-moves) that
    // are constrained (part of a cycle with hitting number 1). mds is assumed
    // to be sorted.
    template<typename C>
    void fill_constrained(const C& mds) {
        for(auto& c : constrained)
            c.reset();

        // Not necessary anymore? Deal with more complicated case in imoves().
        // First mark all the left companions of the homopolymers as
        // constrained.
        // for(mer_t b = 0; b < mer_op_t::alpha; ++b) {
        //     const auto homopoly = mer_op_t::homopolymer(b);
        //     for(mer_t lb = 0; lb < mer_op_t::alpha; ++lb) {
        //         constrained[mer_op_t::lc(homopoly, lb)] = true;
        //     }
        // }
        visit_all(mds);
    }

    void visit_all(const std::vector<mer_t>& mds) {
        for(const auto mer : mds) {
            std::fill(visited.begin(), visited.end(), nil);
            visit(mer, mer, mds, false);
        }
    }

    void visit_all(const std::vector<tristate_t>& bmds) {
         // mds.clear();
        for(mer_t mer = 0; mer < mer_op_t::nb_mers; ++mer) {
            if(includes(bmds, mer)) {
                // mds.push_back(mer);
                std::fill(visited.begin(), visited.end(), nil);
                visit(mer, mer, bmds, false);
            }
        }
    }

    inline static bool includes(const std::vector<mer_t>& mds, const mer_t m) {
        return std::binary_search(mds.cbegin(), mds.cend(), m);
    }

    inline static bool includes(const std::vector<tristate_t>& mds, const mer_t m) {
        return mds[m] == yes;
    }

    template<typename C>
    bool visit(mer_t start, mer_t mer, const C& mds, bool used_fmove) {
        // std::cout << "visit " << start << ' ' << mer << ' ' << used_fmove << std::endl;
        if(visited[mer] != nil && used_fmove) return visited[mer] == yes;

        bool res = false;
        const auto nmer = mer_op_t::nmer(mer);
        const auto cycling = mer_op_t::rb(nmer); // Using that base is staying on PCR

        // std::cout << "\tinfo " << nmer << ' ' << cycling << std::endl;

        for(mer_t b = 0; b < mer_op_t::alpha; ++b) {
            // const auto b = (cycling + i) % mer_op_t::alpha;
            bool ufm = b != cycling; // Traversing a FM
            const mer_t m = mer_op_t::rc(nmer, b);
            // std::cout << "\tloop " << m << ' ' << b << ' ' << ufm << std::endl;
            if(m == start) {
                if(used_fmove) {
                    if(ufm) constrained[mer_op_t::fmove(mer)].set(mer_op_t::lb(mer), b);
                    res = true;
                }
            } else if(!includes(mds, m)) {
                // std::cout << "\trecurse " << m << ' ' << b << ' ' << ufm << std::endl;
                const bool nres = visit(start, m, mds, used_fmove || ufm);
                res = res || nres;
                if(nres && ufm) constrained[mer_op_t::fmove(mer)].set(mer_op_t::lb(mer), b);
            }
        }

        if(used_fmove) visited[mer] = res ? yes : no;
        return res;
    }

    template<typename C>
    imove_sig_t imoves(const C& mds) {
        imove_sig_t res;
        imoves(mds, res);
        return res;
    }

    // I-moves are set of unconstrained edges
    template<typename C>
    void imoves(const C& mds, imove_sig_t& res) {
        res.clear();
        fill_constrained(mds);

        // For combination of F-moves and mask, find those that are not
        // constrained. If f is the homopolymer f = i^(k-1) with i \in \Sigma,
        // then (f, m) with m_i = 1 and all other bits are 0 is not an I-move.
        // It is no move at all. imilarly, (f, m) with m_i = 0 and all other
        // bits are 1 is not an I-move, it is an (degenarated) F-move.
        mer_t base = 0;
        mask_t test1 = (mask_t)1 << base;
        mask_t test2 = ~test1;
        auto nhomo = mer_op_t::fmove(mer_op_t::homopolymer(base));

        for(mer_t fm = 0; fm < mer_op_t::nb_fmoves; ++fm) {
            for(mask_t mask = 1; mask < imove_t::all; ++mask) {
                if(fm == nhomo &&
                   ( ((mask | test1) == imove_t::all) || ((mask & test2) == 0) ) )
                    continue; // One of the degenerated cases: not an I-move

                // Check if any of the bits set in the mask contradict the
                // constrained cycles.
                mask_t mi = 1;
                bool is_possible = true;
                for(mer_t i = 0; is_possible && i < mer_op_t::alpha; ++i, mi <<= 1) {
                    if((mask & mi) == 0) continue;
                    mer_t mj = 1;
                    for(mer_t j = 0; j < mer_op_t::alpha; ++j, mj <<= 1) {
                        if(i == j || (mask & mj ) != 0) continue;
                        if(constrained[fm].get(i, j)) {
                            is_possible = false;
                            break;
                        }
                    }
                }

                if(is_possible)
                    res.emplace_back(fm, mask);
            }

            if(fm == nhomo) {
                ++base;
                nhomo = mer_op_t::fmove(mer_op_t::homopolymer(base));
                test1 = (mask_t)1 << base;
                test2 = ~test1;
            }
        }
    }
};

#endif // _IMOVES_H_
