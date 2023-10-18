#include <cstdlib>
#include <random>
#include <vector>
#include <deque>
#include <bitset>
#include <set>
#include <map>

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

    uint32_t total_imoves = 0;
    mer_t total_moves = 0;
    std::bitset<mer_ops::nb_fmoves> done_moves;
    std::map<std::pair<mer_t,uint8_t>, uint32_t> imoves_use;

    // Initialize at random
    // std::cout << "initial: " << mds << std::endl;

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
            // std::cout << "f-move #" << (size_t)total_moves << ' ' << (size_t)fm << ": " << mds << std::endl;
        }
        // Found an MDS?
        if(total_moves == mer_ops::nb_fmoves) {
            assert2(!mds.fms.empty(), "Found MDS with empty possible F-moves");
            done = found_mds = true;
            break;
        }

        assert2(mds.fms.empty(), "PCR set should have exhausted F-moves");

        // Find I-moves while doing RF-moves
        bool found_imove = false;
        while(!found_imove) {
            while(!mds.rfms.empty() && mds.ims.empty()) {
                const auto fm = *random_set_elt(mds.rfms, rng);
                mds.do_rfmove(fm);
                // std::cout << "rf-move " << (size_t)fm << ": " << mds << std::endl;
            }

            if(mds.ims.empty()) {
                done = true;
                break; // Not found an MDS
            }

            auto ims_it = random_set_elt(mds.ims, rng);
            const auto ims = *ims_it;
            const auto mask = mds.imove_mask(ims);
            const auto ims_uses = ++imoves_use[std::make_pair(ims, mask)];
            found_imove = (ims_uses == 1) || (ims_uses * rng() < 1.0);
            if(found_imove) {
                mds.do_imove(ims);
                ++total_imoves;
                // std::cout << "i-move #" << total_imoves << ' ' << (size_t)ims << ':' << (size_t)mask << ": " << mds <<  std::endl;
                if(total_imoves >= args.max_arg)
                    done = true;
            } else {
                // std::cout << "skip " << (size_t)ims << ':' << (size_t)mask << ' ' << ims_uses << std::endl;
                mds.ims.erase(ims_it);
            }
        }
        // mds after I-move is new starting point. Start over.
    }

    // if(found_mds)
    std::cout << (found_mds ? "mds" : "not") << " : " << mds << std::endl;

    return EXIT_SUCCESS;
}
