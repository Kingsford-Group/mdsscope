#ifndef _IMOVES_H_
#define _IMOVES_H_

#include <vector>
#include <algorithm>

#include "common.hpp"
#include "mer_op.hpp"
#include "mds_op.hpp"

// Figure out the possible I-moves by doing DFSs. Caches data: not multi-thread
// safe (should have 1 object per-thread).
template<typename mer_op_type>
struct imoves_type {
    typedef mer_op_type mer_op_t;
    std::vector<bool> constrained;
    std::vector<tristate_t> visited;
    std::vector<mer_type> mds;
    std::vector<bool> done;

    imoves_type()
        : constrained(mer_op_t::nb_mers)
        , visited(mer_op_t::nb_mers)
        , mds(mer_op_t::nb_mers)
        , done(mer_op_t::nb_fmoves)
    { }

    // Fill the constrained vector marking the edges (potential i-moves) that
    // are constrained (part of a cycle with hitting number 1). mds is assumed
    // to be sorted.
    template<typename C>
    void fill_constrained(const C& mds) {
        std::fill(constrained.begin(), constrained.end(), false);

        // First mark all the left companions of the homopolymers as
        // constrained.
        for(mer_type b = 0; b < mer_op_t::alpha; ++b) {
            const auto homopoly = mer_op_t::homopolymer(b);
            for(mer_type lb = 0; lb < mer_op_t::alpha; ++lb) {
                constrained[mer_op_t::lc(homopoly, lb)] = true;
            }
        }
        visit_all(mds);
    }

    void visit_all(const std::vector<mer_type>& mds) {
        for(const auto mer : mds) {
            std::fill(visited.begin(), visited.end(), nil);
            visit(mer, mer, mds, false);
        }
    }

    void visit_all(const std::vector<tristate_t>& bmds) {
        mds.clear();
        for(mer_type mer = 0; mer < mer_op_t::nb_mers; ++mer) {
            if(includes(bmds, mer)) {
                mds.push_back(mer);
                std::fill(visited.begin(), visited.end(), nil);
                visit(mer, mer, bmds, false);
            }
        }
    }

    inline static bool includes(const std::vector<mer_type>& mds, const mer_type m) {
        return std::binary_search(mds.cbegin(), mds.cend(), m);
    }

    inline static bool includes(const std::vector<tristate_t>& mds, const mer_type m) {
        return mds[m] == yes;
    }

    template<typename C>
    bool visit(mer_type start, mer_type mer, const C& mds, bool used_fmove) {
        if(visited[mer] != nil && used_fmove) return visited[mer] == yes;

        bool res = false;
        const auto nmer = mer_op_t::nmer(mer);
        const auto cycling = mer_op_t::rb(nmer); // Using that base is staying on PCR

        for(mer_type b = 0; b < mer_op_t::alpha; ++b) {
            // const auto b = (cycling + i) % mer_op_t::alpha;
            bool ufm = b != cycling; // Traversing a FM
            const mer_type m = mer_op_t::rc(nmer, b);
            if(m == start) {
                if(used_fmove) {
                    if(ufm) constrained[mer] = true;
                    res = true;
                }
            } else if(!includes(mds, m)) {
                const bool nres = visit(start, m, mds, used_fmove || ufm);
                res = res || nres;
                if(nres && ufm) constrained[mer] = true;
            }
        }

        if(used_fmove) visited[mer] = res ? yes : no;
        return res;
    }

    template<typename C>
    std::vector<mer_type> imoves(const C& mds) {
        std::vector<mer_type> res;
        imoves(mds, res);
        return res;
    }

    // I-moves are unconstrained edges
    template<typename C>
    void imoves(const C& mds, std::vector<mer_type>& res) {
        res.reserve(mds.size());
        fill_constrained(mds);
        for(mer_type i = 0; i < constrained.size(); ++i) {
            if(!constrained[i])
                res.push_back(i);
        }
    }
};

#endif // _IMOVES_H_