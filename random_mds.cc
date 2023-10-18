#include <cstdlib>
#include <iostream>
#include <queue>
#include <stack>
#include <vector>
#include <bitset>
#include <algorithm>

#include <random_mds.hpp>
#include "random_seed.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
#include "mds_op.hpp"
#include "dbg.hpp"
#include "backtrace.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef mds_op_type<mer_ops> mds_ops;

enum fm_type {
  FMNo, // Not an F-move, missing LC
  FMInc, // Incomplete F-move: at least one LC, all others are nil
  FMComp // Complete F-move: all LC are yes
};
fm_type classify_fm(mer_t fm, const std::vector<tristate_t>& bmds) {
    mer_t b = 0;
    bool potential = false;
    for( ; b < mer_ops::alpha; ++b) {
        const auto v = bmds[mer_ops::lc(fm, b)];
        if(v == no) return FMNo;
        potential = potential || (v == nil);
    }
    return potential ? FMInc : FMComp;
}

// Do F-move fm and update bmds accordingly. All the left companion of fm should
// be set to yes or nil in bmds.
//
// New possible F-moves are queued in possible_fms. Any potential F-move (mix of
// yes and nil left-companions) are added to potential_fms. Newly discovered
// pcrs (a nil node becomes yes) are added to new_pcr.
void do_fmove(mer_t fm, std::vector<tristate_t>& bmds,
              std::queue<mer_t>& possible_fms, std::vector<mer_t>& potential_fms, std::vector<mer_t>& new_pcrs) {
    //assert2(mer_op_t::nb_mers == std::accumulate(bmds.cbegin(), bmds.cend(), 0, [](mer_t a, tristate_t x) { return a + (x == yes); }), "Too few mers in bmds");
    // 1) Do the actual F-move and updates PCRs as needed
    for(mer_t b = 0; b < mer_ops::alpha; ++b) {
        const auto m = mer_ops::lc(fm, b);
        auto nm = mer_ops::nmer(m);
        assert2(bmds[m] != no, "Left companion not present " << fm << ' ' << b);
        const bool ftime = bmds[m] == nil; // First time touching this PCR. Mark all k-mers
        bmds[m] = no;
        assert2(bmds[nm] != yes, "Right companion present " << fm << ' ' << b);
        bmds[nm] = yes;

        if(ftime) {
            new_pcrs.push_back(m);
            for(nm = mer_ops::nmer(nm); nm != m; nm = mer_ops::nmer(nm)) {
                assert2(bmds[nm] == nil, "Mix of nil and non-nil mers in PCR: " << m << ' ' << nm);
                bmds[nm] = no;
            }
        }
    }

    // 2) Find potential F-moves. Must be done after the F-move has been
    // completed
    for(mer_t b = 0; b < mer_ops::alpha; ++b) {
        const auto nm = mer_ops::nmer(fm, b);
        const auto nfm = mer_ops::fmove(nm);
        switch(classify_fm(nfm, bmds)) {
        case FMInc: potential_fms.push_back(nfm); break;
        case FMComp: possible_fms.push(nfm); break;
        default: break;
        }

    }
    // assert2(mer_op_t::nb_mers == std::accumulate(bmds.cbegin(), bmds.cend(), 0, [](mer_t a, tristate_t x) { return a + (x == yes); }), "Too few mers in bmds");
}

// Undo an F-move fm and update bmds accordingly. All the right companion of fm
// should be set to yes in bmds. Calling undo_fmove after do_fmove should return
// bmds in its pre do_fmove state.
//
// new_pcr should be as set by do_fmove. All k-mers in these pcrs are set to nil
// again.
void undo_fmove(mer_t fm, std::vector<tristate_t>& bmds, const std::vector<mer_t>& new_pcrs) {
    for(mer_t b = 0; b < mer_ops::alpha; ++b) {
        const auto m = mer_ops::nmer(fm, b);
        auto pm = mer_ops::pmer(m);
        assert2(bmds[m] != no, "Right companion not present " << fm << ' ' << b);
        bmds[m] = no;
        assert2(bmds[pm] != yes, "Left companion present " << fm << ' ' << b);
        bmds[pm] = yes;
        if(std::find(new_pcrs.begin(), new_pcrs.end(), pm) != new_pcrs.end()) {
            // Reset PCR back to nil
            for( ; pm != m; pm = mer_ops::pmer(pm))
                bmds[pm] = nil;
            bmds[m] = nil;
        }
    }
}

// Efficient DFS. Only start from subset of nodes.
template<typename mer_op_type>
struct partial_dfs_type {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_t::mer_t mer_t;

    std::bitset<mer_op_t::nb_mers> visiting;
    std::bitset<mer_op_t::nb_mers> visited;
    std::stack<std::pair<mer_t, mer_t>> stack; // first: mer, second: base

    bool has_cycle(const std::vector<tristate_t>& bmds, const std::vector<mer_t>& new_pcrs) {
        visiting.reset();
        visited.reset();
        // for(const auto pcr : new_pcrs) {
        //     // Start from node after the selected one in the PCR
        //     auto start = pcr;
        //     if(bmds[start] != yes) {
        //         for(start = mer_op_t::nmer(start); start != pcr && bmds[start] != yes; start = mer_op_t::nmer(start)) ;
        //     }
        //     start = mer_op_t::nmer(start);
        //     assert2(start == pcr || bmds[start] != yes, "Didn't find proper start: start shouldn't be selected" << pcr << ' ' << start);
        //     assert2(bmds[mer_op_t::pmer(start)] == yes, "Didn't find proper start: node before start should be selected " << pcr << ' ' << start);
        //     if(bmds[start] != yes && visit(bmds, start))
        //         return true;
        // }
        for(mer_t start = 0; start < mer_op_t::nb_mers; ++start) {
            if(bmds[start] == no && visit(bmds, start))
                return true;
        }
        return false;
    }

    bool visit(const std::vector<tristate_t>& bmds, mer_t start) {
        stack.emplace(start, 0);
        visiting[start] = true;

        mer_t node, b;
        while(!stack.empty()) {
            std::tie(node, b) = stack.top();
            stack.pop();
            if(b >= mer_op_t::alpha) {
                visiting.set(node, false);
                visited.set(node);
                continue;
            }
            stack.emplace(node, b + 1);

            auto nnode = mer_op_t::nmer(node, b);
            if(visiting[nnode])
                return true; // Found back edge
            if(!visited[nnode] && bmds[nnode] == no) {
                stack.emplace(nnode, 0);
                visiting.set(nnode);
            }
        }

        return false;
    }
};
typedef partial_dfs_type<mer_ops> partial_dfs_t;

int main(int argc, char* argv[]) {
    show_backtrace();
    random_mds args(argc, argv);

    auto prg = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
                                           args.iseed_given ? args.iseed_arg : nullptr);


    std::vector<tristate_t> mds(mer_ops::nb_mers, nil);
    std::vector<tristate_t> saved_mds;
    std::vector<bool> done_fms(mer_ops::nb_fmoves, false);
    mer_t done_total = 0;
    partial_dfs_t dfs;

    // Add homopolymers to mds
    for(mer_t b = 0; b < mer_ops::alpha; ++b)
        mds[mer_ops::homopolymer(b)] = yes;

    std::uniform_int_distribution<mer_t> fmove_rng(0, mer_ops::nb_fmoves);
    const mer_t initial_fm = fmove_rng(prg);
    std::cout << "Initial F-move " << (size_t)initial_fm << std::endl;

    std::queue<mer_t> possible_fms, saved_possible_fms;
    std::vector<mer_t> potential_fms, saved_potential_fms;
    std::vector<mer_t> new_pcrs;
    possible_fms.push(initial_fm);

    while(done_total < mer_ops::nb_fmoves && !(possible_fms.empty() && potential_fms.empty())) {
        // 1) Apply all possible F-moves
        while(done_total < mer_ops::nb_fmoves && !possible_fms.empty()) {
            const auto fm = possible_fms.front();
            possible_fms.pop();

            new_pcrs.clear();
            std::cout << "Apply F-move " << (size_t)fm << std::endl;;
            do_fmove(fm, mds, possible_fms, potential_fms, new_pcrs);
            if(!done_fms[fm]) {
                done_fms[fm] = true;
                ++done_total;
            }
        }

        // 2) Find a valid potential F-move and do it. Repeat until other possible fmoves
        while(done_total < mer_ops::nb_fmoves && !potential_fms.empty()) {
            const auto fm = potential_fms.back();
            potential_fms.pop_back();

            if(classify_fm(fm, mds) == FMInc) {
                new_pcrs.clear();
                std::cout << "Apply potential F-move " << (size_t)fm << std::endl;
                saved_mds = mds;
                saved_possible_fms = possible_fms;
                saved_potential_fms = potential_fms;
                do_fmove(fm, mds, possible_fms, potential_fms, new_pcrs);
                std::cout << "New PCRs: " << joinT<size_t>(new_pcrs, ',') << std::endl;
                if(dfs.has_cycle(mds, new_pcrs)) {
                    std::cout << "Undo F-move " << (size_t)fm << std::endl;
                    mds = saved_mds;
                    possible_fms = saved_possible_fms;
                    potential_fms = saved_potential_fms;
                    // undo_fmove(fm, mds, new_pcrs);
                    // assert2(saved_mds == mds, "Undo didn't work properly");
                    continue;
                }
                done_fms[fm] = true;
                ++done_total;
                break;
            }
        }
    }

    int ret = EXIT_SUCCESS;
    if(done_total < mer_ops::nb_fmoves) {
        std::cerr << "Stock in non-decycling PCR set " << (size_t)done_total << ' ' << (size_t)mer_ops::nb_fmoves << std::endl;
        ret = EXIT_FAILURE;
    }

    std::cout << mds << '\n';

    return ret;
}
