#include <cstdlib>
#include <random>
#include <stack>
#include <vector>
#include <deque>
#include <bitset>
#include <set>
#include <map>
#include <algorithm>

#include "random_pcr.hpp"

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
    NMove, // Nothing
    FMove,
    RFMove,
    IMove,
};

template<size_t S>
fm_type classify_fm(mer_t fm, const std::bitset<S>& bmds) {
    unsigned int lc = 0, rc = 0;
    const mer_t rfm = fm * mer_ops::alpha;
    for(mer_t b = 0; b < mer_ops::alpha; ++b) {
        if(bmds.test(mer_ops::lc(fm, b)))
            ++lc;
        if(bmds.test(mer_ops::rc(rfm, b)))
            ++rc;
    }
    // rc + lc can be > alpha as homopolymers count for both lc and rc
    const bool is_homo_fm = mer_ops::is_homopolymer_fm(fm);
    assert2(lc + rc <= mer_ops::alpha + (is_homo_fm ? 1 : 0), "Classify fm found more lc + rc than alpha");
    // std::cout << "classify " << (size_t)fm << ' ' << (size_t)rfm << ": " << joinT<size_t>(bmds, ',') << " | " << lc << ' ' << rc << '\n';

    if(lc == mer_ops::alpha) return FMove;
    if(rc == mer_ops::alpha) return RFMove;
    if(lc > 0 && rc > 0 && (rc + lc == mer_ops::alpha + (is_homo_fm ? 1 : 0))) return IMove;
    return NMove;
}

struct mds_fms {
    std::bitset<mer_ops::nb_mers> bmds;
    std::set<mer_t> fms, rfms, ims;

    template<typename R>
    mds_fms(R& rng) {
        for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
            unsigned int pcr_size = 1;
            for(mer_t nm = mer_ops::nmer(m); nm != m && pcr_size > 0; nm = mer_ops::nmer(nm)) {
                if(nm < m)
                    pcr_size = 0;
                else
                    ++pcr_size;
            }
            if(pcr_size == 0) continue; // Only look at PCRs once

            const auto choose = std::uniform_int_distribution<mer_t>(0, pcr_size-1)(rng);
            mer_t nm = m;
            for(mer_t i = 0; i < choose; ++i)
                nm = mer_ops::nmer(nm);
            bmds.set(nm);
        }

        for(mer_t fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
            switch(classify_fm(fm, bmds)) {
            case FMove: fms.insert(fm); break;
            case RFMove: rfms.insert(fm); break;
            case IMove: ims.insert(fm); break;
            default: break;
            }
        }
    }

    void do_fmove(const mer_t fm) {
        assert2(fms.find(fm) != fms.end(), "Not a possible F-move " << fm);
        assert2(rfms.find(fm) == rfms.end(), "Shouldn't be a possible RF-move " << fm);
        assert2(ims.find(fm) == ims.end(), "Shouldn't be a possible I-move " << fm);

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto m = mer_ops::lc(fm, b);
            const auto nm = mer_ops::nmer(m);
            assert2(bmds.test(m), "Left companion not present " << fm << ' ' << b);
            bmds.flip(m);
            assert2(!bmds.test(m), "Right companion present " << fm << ' ' << b);
            bmds.flip(nm);
        }

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto nfm = mer_ops::fmove(mer_ops::nmer(fm, b));
            switch(classify_fm(nfm, bmds)) {
            case FMove:
                fms.insert(nfm);
                ims.erase(nfm);
                break;
            case IMove:
                ims.insert(nfm);
                break;
            case RFMove:
                assert2(mer_ops::is_homopolymer_fm(nfm), "Doing an F-move created a rogue RF-move " << fm << ' ' << nfm);
                break;
            default: break;
            }

            const auto pfm = mer_ops::fmove(mer_ops::pmer(mer_ops::lc(fm, b)));
            rfms.erase(pfm);
            ims.erase(pfm);
        }

        fms.erase(fm);
        rfms.insert(fm);
    }

    void do_rfmove(const mer_t fm) {
        assert2(rfms.find(fm) != rfms.end(), "Not a possible RF-move " << fm);
        assert2(fms.find(fm) == fms.end(), "Shouldn't be a possible F-move " << fm);
        assert2(ims.find(fm) == ims.end(), "Shouldn't be a possible I-move " << fm);

        const mer_t rfm = fm * mer_ops::alpha;
        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto m = mer_ops::rc(rfm, b);
            const auto pm = mer_ops::pmer(m);
            assert2(bmds.test(m), "Right companion not present " << rfm << ' ' << b);
            bmds.flip(m);
            assert2(!bmds.test(pm), "Left companion present " << rfm << ' ' << b);
            bmds.flip(pm);
        }

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto pfm = mer_ops::fmove(mer_ops::pmer(mer_ops::lc(fm, b)));
            switch(classify_fm(pfm, bmds)) {
            case RFMove:
                rfms.insert(pfm);
                ims.erase(pfm);
                break;
            case IMove:
                ims.insert(pfm);
                break;
            case FMove:
                assert2(mer_ops::is_homopolymer_fm(pfm), "Doing an RF-move created a rogue F-move " << (size_t)fm << ' ' << (size_t)b << ' ' << (size_t)pfm);
                break;
            default: break;
            }

            const auto nfm = mer_ops::fmove(mer_ops::nmer(fm, b));
            fms.erase(nfm);
            ims.erase(nfm);
        }

        rfms.erase(fm);
        fms.insert(fm);
    }

    void do_imove(const mer_t fm) {
        assert2(ims.find(fm) != ims.end(), "Not a possible I-move " << fm);
        assert2(fms.find(fm) == fms.end(), "Shouldn't a possible F-move " << fm);
        assert2(rfms.find(fm) == rfms.end(), "Shouldn't be a possible RF-move " << fm);

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto m = mer_ops::lc(fm, b);
            const auto nm = mer_ops::nmer(m);
            if(m != nm) {
                assert2(bmds.test(m) ^ bmds.test(nm), "Neither left or right companion present " << (size_t)fm << ' ' << b);
                bmds.reset(m);
                bmds.set(nm);
            }
        }

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto nfm = mer_ops::fmove(mer_ops::nmer(fm, b));
            switch(classify_fm(nfm, bmds)) {
            case FMove:
                fms.insert(nfm);
                ims.erase(nfm);
                break;
            case IMove:
                ims.insert(nfm);
                break;
            case RFMove:
                assert2(mer_ops::is_homopolymer_fm(nfm), "Doing an I-move created a rogue RF-move " << (size_t)nfm << ' ' << b);
                break;
            default: break;
            }

            const auto pfm = mer_ops::fmove(mer_ops::pmer(mer_ops::lc(fm, b)));
            rfms.erase(pfm);
            ims.erase(pfm);
        }

        ims.erase(fm);
        rfms.insert(fm);
    }

    uint8_t imove_mask(const mer_t im) {
        uint8_t res = 0;
        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            if(bmds.test(mer_ops::lc(im, b)))
                res |= (uint8_t)1 << b;
        }
        return res;
    }

    void move_mer(const mer_t m) {
        assert2(bmds.test(m), "Moving mer valid only for selected mers");
        assert2(!mer_ops::is_homopolymer(m), "Move mer not well defined for homopolymers");

        const auto nm = mer_ops::nmer(m);
        const auto fm = mer_ops::fmove(m);
        const auto pfm = mer_ops::fmove(mer_ops::pmer(m));
        const auto nfm = mer_ops::fmove(nm);

        fms.erase(fm);
        rfms.erase(pfm);
        ims.erase(pfm);

        bmds.reset(m);
        bmds.set(nm);

        switch(classify_fm(fm, bmds)) {
        case RFMove:
            rfms.insert(fm);
            ims.erase(fm);
            break;
        case IMove:
            ims.insert(fm);
            break;
        case FMove:
            assert2(false, "move mer created rogue F-move " << (size_t)m);
            break;
        default: break;
        }

        switch(classify_fm(nfm, bmds)) {
        case FMove:
            fms.insert(nfm);
            ims.erase(nfm);
            break;
        case IMove:
            ims.insert(nfm);
            break;
        case RFMove:
            assert2(false, "move mer created rogue RF-move " << (size_t)m);
            break;
        default: break;

        }
    }
};

template<size_t S>
struct bitset_iterator {
    typedef ssize_t difference_type;
    size_t i;
    const std::bitset<S> *s;

    bitset_iterator(const std::bitset<S>& set) : i(0), s(&set) {
        for(i = 0; i < S && !s->test(i); ++i) ;
    }
    bitset_iterator() : i(S), s(nullptr) {}

    size_t operator*() const { return i; }
    bool operator==(const bitset_iterator& rhs) const { return i == rhs.i; }
    bool operator!=(const bitset_iterator& rhs) const { return i != rhs.i; }
    bitset_iterator& operator++() {
        for(++i; i < S && !s->test(i); ++i) ;
        return *this;
    }
    bitset_iterator operator++(int) {
        bitset_iterator ret(*this);
        ++*this;
        return ret;
    }
};

namespace std {
template<size_t S>
bitset_iterator<S> begin(const std::bitset<S>& set) { return bitset_iterator<S>(set); }
template<size_t S>
bitset_iterator<S> end(const std::bitset<S>& set) { return bitset_iterator<S>(); }
}

std::ostream& operator<<(std::ostream& os, const mds_fms& mds) {
    return os << '{'
              << joinT<size_t>(mds.bmds, ',') << " F "
              << joinT<size_t>(mds.fms, ',') << " R "
              << joinT<size_t>(mds.rfms, ',') << " I "
              << joinT<size_t>(mds.ims, ',')
              << '}';
}

std::ostream& operator<<(std::ostream& os, const std::pair<const mer_t, uint8_t>& x) {
    return os << (size_t)x.first << ':' << (size_t)x.second;
}

std::ostream& operator<<(std::ostream& os, const std::pair<mer_t, mer_t>& x) {
    return os << (size_t)x.first << ':' << (size_t)x.second;
}


template<typename T, typename R>
typename std::set<T>::iterator random_set_elt(std::set<T>& s, R& rng) {
    auto ret = s.begin();
    const auto size = s.size();
    if(size > 1) {
        const auto skip = std::uniform_int_distribution<decltype(size)>(0, size - 1)(rng);
        std::advance(ret, skip);
    }
    return ret;
}

template<typename R>
mer_t find_special_cycle(const std::bitset<mer_ops::nb_mers>& mds, std::vector<std::pair<mer_t, mer_t>>& scycle, R& rng) {
    scycle.clear();

    // Start with a random mer that is not a homopolymer
    std::uniform_int_distribution<mer_t> mer_rng(0, mer_ops::nb_mers-1);
    mer_t mer = mer_rng(rng);
    while(mer_ops::is_homopolymer(mer))
        mer = mer_rng(rng);

    mer_t lcs[mer_ops::alpha]; // List of non-selected left-companions

    // Traverse the PCR backward (following pmer()) until the selected mer of
    // the PCR. If this mer is already in scycle, we found the cycle and we are
    // done. Otherwise, add to scycle and choose at random one of the
    // non-selected left companion (which must exists as there are no possible
    // F-move).
    while(true) {
        for( ; !mds.test(mer); mer = mer_ops::pmer(mer)) ;

        auto it = std::find_if(scycle.begin(), scycle.end(), [&](const std::pair<mer_t, mer_t>& m) { return m.first == mer; });
        if(it != scycle.end())
            return std::distance(scycle.begin(), it); // Found a cycle. Returned offset is start of cycle

        unsigned nb_lcs = 0; // Number of unselected lc of mer found
        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            mer_t lcm = mer_ops::lc(mer, b);
            if(!mds.test(lcm))
                lcs[nb_lcs++] = lcm;
        }
        assert2(nb_lcs > 0, "No unselected left-companion in set with no possible F-move");

        const auto nexti = std::uniform_int_distribution<unsigned>(0, nb_lcs-1)(rng);
        scycle.emplace_back(mer, lcs[nexti]);
        mer = lcs[nexti];
    }

    // Should never get here
    return mer_ops::nb_mers;
}

// DFS from source to source. Used to check for the existence of weight 1 cycles
// that are not PCRs.
// struct partial_dfs_type {
//     std::bitset<mer_ops::nb_mers> visiting;
//     std::bitset<mer_ops::nb_mers> visited;
//     std::stack<std::pair<mer_t, mer_t>> stack; // first: mer, second: base

//     bool has_cycle(const std::bitset<mer_ops::nb_mers>& bmds, const mer_t start) {
//         visiting.reset();
//         visited.reset();


//         for(mer_t start = 0; start < mer_ops::nb_mers; ++start) {
//             if(bmds[start] == no && visit(bmds, start))
//                 return true;
//         }
//         return false;
//     }

//     bool visit(const std::vector<tristate_t>& bmds, mer_t start) {
//         stack.emplace(start, 0);
//         visiting[start] = true;

//         mer_t node, b;
//         while(!stack.empty()) {
//             std::tie(node, b) = stack.top();
//             stack.pop();
//             if(b >= mer_ops::alpha) {
//                 visiting.set(node, false);
//                 visited.set(node);
//                 continue;
//             }
//             stack.emplace(node, b + 1);

//             auto nnode = mer_ops::nmer(node, b);
//             if(visiting[nnode])
//                 return true; // Found back edge
//             if(!visited[nnode] && bmds[nnode] == no) {
//                 stack.emplace(nnode, 0);
//                 visiting.set(nnode);
//             }
//         }

//         return false;
//     }
// };


// Algo:
//
// Start from a random PCR M.
// Then, Main Loop:
//   Do all possible F-moves.
//   If can do all \sigma^(k-1) F-moves, output M, done. (it is an MDS)
//   If not, do all possible RF-moves looking for a possible I-moves
//   If find I-move: do I-move, update M accordingly, repeat main loop.
//   If exhaust RF-moves: fail

int main(int argc, char* argv[]) {
    random_pcr args(argc, argv);

    auto rng = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
                                           args.iseed_given ? args.iseed_arg : nullptr);

    struct mds_fms mds(rng);

    mer_t total_moves = 0;
    std::bitset<mer_ops::nb_fmoves> done_moves;

    std::vector<std::pair<mer_t, mer_t>> scycle; // Special hitting number 0 cycle
    uint32_t mer_moves = 0;

    // Initialize at random
    std::cout << "initial: " << mds << std::endl;

    auto found_mds = false;
    auto done = false;
    while(!done) {
        // Do all possible F-moves
        total_moves = 0;
        done_moves.reset();
        while(!mds.fms.empty() && total_moves < mer_ops::nb_fmoves) {
            const auto fm = *mds.fms.begin();
            mds.do_fmove(fm);
            if(!done_moves.test(fm)) {
                done_moves.set(fm);
                ++total_moves;
            }
            std::cout << "f-move #" << (size_t)total_moves << ' ' << (size_t)fm << ": " << mds << std::endl;
        }
        // Found an MDS?
        if(total_moves == mer_ops::nb_fmoves) {
            assert2(!mds.fms.empty(), "Found MDS with empty possible F-moves");
            done = found_mds = true;
            break;
        }

        assert2(mds.fms.empty(), "PCR set should have exhausted F-moves");

        // Find a special hitting number zero cycle and fix it. This cycle is
        // guaranteed to exists because there are no possible F-move. Then pick
        // a random mer on the cycle and move it into the cycle.
        const auto offset = find_special_cycle(mds.bmds, scycle, rng);
        std::cout << "scycle " << (size_t)offset << ' ' << join(scycle, ',') << std::endl;
        assert2(!scycle.empty(), "Special cycle is empty");
        assert2(offset < scycle.size(), "Start offset is larger than scycle.size");
        const auto movei = std::uniform_int_distribution<size_t>(offset, scycle.size() - 1)(rng);
        assert2(mds.bmds.test(scycle[movei].first), "Mer to move is not selected");
        mds.move_mer(scycle[movei].first);
        ++mer_moves;
        if(mer_moves >= args.max_arg)
            done = true;
        std::cout << "m-move #" << (size_t)mer_moves << ' ' << (size_t)scycle[movei].first << ": " << mds << std::endl;

        // mds after mer move is new starting point. Start over.
    }

    // if(found_mds)
    std::cout << (found_mds ? "mds" : "not") << " : " << mds << std::endl;

    return EXIT_SUCCESS;
}
