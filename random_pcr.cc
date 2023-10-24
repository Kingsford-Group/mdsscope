#include <cstdlib>
#include <random>
#include <stack>
#include <vector>
#include <deque>
#include <bitset>
#include <set>
#include <map>
#include <algorithm>
#include <queue>

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
typedef std::bitset<mer_ops::nb_mers> bmds_t;
typedef std::bitset<mer_ops::nb_fmoves> bfms_t;

enum fm_type {
    NMove, // Nothing
    FMove,
    RFMove,
    IMove,
};

fm_type classify_fm(mer_t fm, const bmds_t& bmds) {
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
    bmds_t bmds;
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

    // Move one mer forward. Not necessarily a signature preserving operation.
    // Needs to be checked independently.
    void fmove_mer(const mer_t m) {
        assert2(bmds.test(m), "F-Moving mer valid only for selected mers");
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
            assert2(false, "fmove mer created rogue F-move " << (size_t)m);
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

    void rmove_mer(const mer_t m) {
        assert2(bmds.test(m), "R-Moving mer valid only for selected mers");
        assert2(!mer_ops::is_homopolymer(m), "R-Move mer not well defined for homopolymer");

        const auto pm = mer_ops::pmer(m);
        const auto fm = mer_ops::fmove(m);
        const auto pfm = mer_ops::fmove(pm);
        const auto ppfm = mer_ops::fmove(mer_ops::pmer(pm));

        fms.erase(fm);
        rfms.erase(pfm);
        ims.erase(fm);

        bmds.reset(m);
        bmds.set(pm);

        switch(classify_fm(pfm, bmds)) {
        case FMove:
            fms.insert(pfm);
            ims.erase(pfm);
            break;
        case IMove:
            ims.insert(pfm);
            break;
        case RFMove:
            assert2(false, "rmove mer created rogue RF-move " << (size_t)m);
            break;
        default: break;
        }

        switch(classify_fm(ppfm, bmds)) {
        case RFMove:
            rfms.insert(ppfm);
            ims.erase(ppfm);
            break;
        case IMove:
            ims.insert(ppfm);
            break;
        case FMove:
            assert2(false, "rmove mer create rogue F-move " << (size_t)m);
            break;
        default: break;
        }
    }

    template<typename PF, typename R>
    mer_t check_imoves(bfms_t& checked_ims, PF& path_finder, R& rng) {
        std::vector<mer_t> shuffled_ims(ims.begin(), ims.end());
        std::shuffle(shuffled_ims.begin(), shuffled_ims.end(), rng);

        for(auto im : shuffled_ims) {
            if(checked_ims.test(im))
                continue; // Already checked, F-moves don't change hitting number so no need to check again
            checked_ims.set(im);
            auto loop_len = path_finder.has_path_lc(bmds, im);
            std::cout << "im loop_len " << (size_t)im << ' ' << (size_t)loop_len << std::endl;
            if(loop_len > 0) { // cycle through I-move. Close it
                do_imove(im);
                return im;
            }
        }

        return mer_ops::nb_fmoves; // sentinel value -> no I-move done
    }
};

struct bitset_iterator {
    typedef ssize_t difference_type;
    size_t i;
    const bmds_t *s;

    bitset_iterator(const bmds_t& set) : i(0), s(&set) {
        for(i = 0; i < mer_ops::nb_mers && !s->test(i); ++i) ;
    }
    bitset_iterator() : i(mer_ops::nb_mers), s(nullptr) {}

    size_t operator*() const { return i; }
    bool operator==(const bitset_iterator& rhs) const { return i == rhs.i; }
    bool operator!=(const bitset_iterator& rhs) const { return i != rhs.i; }
    bitset_iterator& operator++() {
        for(++i; i < mer_ops::nb_mers && !s->test(i); ++i) ;
        return *this;
    }
    bitset_iterator operator++(int) {
        bitset_iterator ret(*this);
        ++*this;
        return ret;
    }
};

namespace std {
bitset_iterator begin(const bmds_t& set) { return bitset_iterator(set); }
bitset_iterator end(const bmds_t& set) { return bitset_iterator(); }
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

// Given a set with no F-move, find a special cycle. The nodes of the set and
// the cycles are added to scycle. Returns the offset into scycle where the
// cycle actually starts. Returns mer_ops::nb_mers if hits mers already checked
// (hence no new cycle found from the starting point mer).
//
// fwd == true means all F-moves have been exhausted. fwd == false, all RF-move
// are exhausted.
template<typename R>
mer_t find_special_cycle(const bmds_t& mds, bool fwd, std::vector<std::pair<mer_t, mer_t>>& scycle, bmds_t& checked_mers, mer_t mer, R& rng) {
    scycle.clear();

    mer_t lcs[mer_ops::alpha]; // List of non-selected left-companions

    // Traverse the PCR backward (following pmer()) until the selected mer of
    // the PCR. If this mer is already in scycle, we found the cycle and we are
    // done. Otherwise, add to scycle and choose at random one of the
    // non-selected left companion (which must exists as there are no possible
    // F-move).
    while(true) {
        for( ; !mds.test(mer); mer = fwd ? mer_ops::pmer(mer) : mer_ops::nmer(mer)) ;

        auto it = std::find_if(scycle.begin(), scycle.end(), [&](const std::pair<mer_t, mer_t>& m) { return m.first == mer; });
        if(it != scycle.end())
            return std::distance(scycle.begin(), it); // Found a cycle. Returned offset is start of cycle

        if(checked_mers.test(mer))
            return mer_ops::nb_mers;

        unsigned nb_cs = 0; // Number of unselected lc of mer found
        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            mer_t cm = fwd ? mer_ops::lc(mer, b) : mer_ops::rc(mer, b);
            if(!mds.test(cm))
                lcs[nb_cs++] = cm;
        }
        assert2(nb_cs > 0, "No unselected left-companion in set with no possible F-move");

        const auto nexti = std::uniform_int_distribution<unsigned>(0, nb_cs-1)(rng);
        scycle.emplace_back(mer, lcs[nexti]);
        checked_mers.set(mer);
        mer = lcs[nexti];
    }

    // Should never get here
    return mer_ops::nb_mers;
}

struct shortest_path {
    bmds_t visited;
    std::deque<std::pair<mer_t, mer_t>> queue; // (mer, distance in PCRs)

    // Assume that the queue has been primed. Do bfs algorithm until predicate
    // say we reach the target mer.
    template<typename P>
    mer_t visit(const bmds_t& bmds, bool fwd, P predicate) {
        mer_t mer, dist;
        while(!queue.empty()) {
            std::tie(mer, dist) = queue.front();
            queue.pop_front();

            for( ; !bmds.test(mer); mer = fwd ? mer_ops::nmer(mer) : mer_ops::pmer(mer)) {
                if(predicate(mer))
                    return dist;
                for(mer_t b = 0; b < mer_ops::alpha; ++b) {
                    if(b == (fwd ? mer_ops::lb(mer) : mer_ops::rb(mer)))
                        continue;
                    mer_t nmer = fwd ? mer_ops::nmer(mer, b) : mer_ops::pmer(mer, b);
                    if(predicate(nmer)) return dist + 1;
                    if(!bmds.test(nmer) && !visited.test(nmer)) {
                        queue.emplace_back(nmer, dist + 1);
                        visited.set(nmer);
                    }
                }
            }
        }

        return 0; // No path found
    }

    // Start from every right-companion of target, except the one on the same
    // PCR, and find a shortest path back to target. Returns the length of the
    // shortest path if any, 0 otherwise.
    mer_t has_path_target(const bmds_t& bmds, bool fwd, const mer_t target) {
        visited.reset();
        queue.clear();

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            if(b == (fwd ? mer_ops::lb(target) : mer_ops::rb(target)))
                continue;
            mer_t mer = fwd ? mer_ops::nmer(target, b) : mer_ops::pmer(target, b);
            if(!bmds.test(mer)) {
                queue.emplace_back(mer, 0);
                visited.set(mer);
            }
        }

        return visit(bmds, fwd, [target](mer_t m) { return m == target; });
    }

    // Start from every right-companion of fm that is not selected. The targets
    // are any non-selected left-companion of fm.
    mer_t has_path_lc(const bmds_t& bmds, const mer_t fm) {
        visited.reset();
        queue.clear();

        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto rc = mer_ops::nmer(fm, b);
            if(!bmds.test(rc)) {
                queue.emplace_back(rc, 0);
                visited.set(rc);
            }
        }
        return visit(bmds, true, [&](mer_t m) { return !bmds.test(m) && mer_ops::fmove(m) == fm; });
    }
};

// Do all possible F-move or RF-moves (depending on the fwd parameter) and
// checkes for I-moves cycles as well. Return 2 values: the number of moves done
// and whether an I-move was done. If and I-move was done, the main loop should
// be restarted. If not and the number of moves done is mer_ops::nb_fmoves, then
// mds is an MDS.
template<typename PF, typename R>
std::pair<mer_t, bool> do_all_possible_moves(mds_fms& mds, bool fwd, PF& path_finder, bfms_t& checked_ims, bfms_t& done_moves, R& rng) {
    bool found_im = false;
    mer_t total_moves = 0;
    done_moves.reset();
    auto& set = fwd ? mds.fms : mds.rfms;

    while(!set.empty() && total_moves < mer_ops::nb_fmoves) {
        // See if doing an I-move would fix a hitting number 0 cycle
        const mer_t im = mds.check_imoves(checked_ims, path_finder, rng);
        if(im < mer_ops::nb_fmoves) {
            std::cout << "i-move " << (size_t)im << ": " << mds << std::endl;
            found_im = true;
            break;
        }

        const auto fm = *random_set_elt(set, rng);
        fwd ? mds.do_fmove(fm) : mds.do_rfmove(fm);
        if(!done_moves.test(fm)) {
            done_moves.set(fm);
            ++total_moves;
        }
        std::cout << (fwd ? "f-move #" : "r-move #") << (size_t)total_moves << ' ' << (size_t)fm << ": " << mds << std::endl;
    }
    if(!found_im) { // Last check for I-moves after last F-moves done (or none)
        const mer_t im = mds.check_imoves(checked_ims, path_finder, rng);
        if(im < mer_ops::nb_fmoves) {
            std::cout << "i-move " << (size_t)im << ": " << mds << std::endl;
            found_im = true;
        }
    }

    return std::make_pair(total_moves, found_im);
}

// Find a special hitting number zero cycle and fix it. This cycle is guaranteed
// to exists when there are no possible F-move or RF-move. Then pick a random
// mer on the cycle that doesnt create another cycle with hitting number 0 when
// moved and move it into the cycle.
//
// If such a mer exists, return true. If not, return false and the mer than when
// moved would create the smallest (in number of PCRs) hitting number 0 cycle.
template<typename PF, typename R>
std::pair<mer_t, bool> do_all_special_cycles(mds_fms& mds, bool fwd, PF& path_finder, bmds_t& checked_mers, R& rng) {
    std::vector<std::pair<mer_t, mer_t>> scycle; // Special hitting number 0 cycle
    auto min_loop = std::numeric_limits<mer_t>::max();
    mer_t min_mer = mer_ops::nb_mers;

    const auto mer_offset = std::uniform_int_distribution<mer_t>(0, mer_ops::nb_mers-1)(rng);
    checked_mers.reset();
    for(mer_t start_i = 0; start_i < mer_ops::nb_mers; ++start_i) {
        const mer_t start = (mer_offset + start_i) % mer_ops::nb_mers;
        if(mer_ops::is_homopolymer(start) || !mds.bmds.test(start) || checked_mers.test(start))
            continue;

        const auto offset = find_special_cycle(mds.bmds, fwd, scycle, checked_mers, start, rng);
        if(offset == mer_ops::nb_mers)
            continue; // Ran into a previously known cycle. Skip
        std::cout << "scycle " << (size_t)offset << ' ' << join(scycle, ',') << std::endl;
        assert2(!scycle.empty(), "Special cycle is empty");
        assert2(offset < scycle.size(), "Start offset is larger than scycle.size");

        auto start_scycle = scycle.begin();
        std::advance(start_scycle, offset);
        std::shuffle(start_scycle, scycle.end(), rng);

        for(auto i = offset; i < scycle.size() && min_loop != 0; ++i) {
            const auto m = scycle[i].first;
            const auto loop_len = path_finder.has_path_target(mds.bmds, fwd, m);
            std::cout << "has path " << (size_t)m << ' ' << (size_t)loop_len << std::endl;
            if(loop_len < min_loop) {
                min_loop = loop_len;
                min_mer = scycle[i].first;
            }
        }

        // Use first found edge with no cycle with hitting number 1
        if(min_loop == 0) {
            fwd ? mds.fmove_mer(min_mer) : mds.rmove_mer(min_mer);
            std::cout << "m-move " << (size_t)min_mer << ": " << mds << std::endl;
            return std::make_pair(min_mer, true);
        }
    }

    return std::make_pair(min_mer, false);
}

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

    mer_t total_moves = 0, min_mer;
    bool done_mer_operation;
    bfms_t done_moves, checked_ims;

    bmds_t checked_mers;
    // uint32_t mer_moves = 0;

    shortest_path path_finder;

    // Initialize at random
    std::cout << "initial: " << mds << std::endl;

    auto found_mds = false;
    auto done = false;
    auto fwd = true; // Alternate between F-moves and RF-moves
    for( ; !done; fwd = !fwd) {
        // Do all possible F-moves / RF-moves
        checked_ims.reset();
        std::tie(total_moves, done_mer_operation) = do_all_possible_moves(mds, fwd, path_finder, checked_ims, done_moves, rng);
        if(done_mer_operation)
            continue;

        if(total_moves == mer_ops::nb_fmoves) {
            assert2(!mds.fms.empty() && !mds.rfms.empty(), "Found MDS with empty possible F-moves or RF-moves");
            done = found_mds = true;
            continue;
        }

        // Do special cycle operations when no F-move / RF-move
        std::tie(min_mer, done_mer_operation) = do_all_special_cycles(mds, fwd, path_finder, checked_mers, rng);
        // if(done_mer_operation)
        //     continue;
    }

    // if(found_mds)
    std::cout << (found_mds ? "mds" : "not") << " : " << mds << std::endl;

    return EXIT_SUCCESS;
}
