#ifndef MDS_OP_H_
#define MDS_OP_H_

#include <set>
#include <algorithm>

#include "dbg.hpp"
#include "mer_op.hpp"
#include "imove_signature.hpp"
#include "common.hpp"

template<typename mer_op_type>
struct mds_op_type {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;
    typedef imove_type<mer_op_type> imove_t;
    typedef typename imove_t::mask_type mask_t;

    std::vector<tristate_t> bmds, nbmds;
    std::vector<bool> done;
    std::vector<mer_t> fmoves, nfmoves, fm_listA;

    mds_op_type()
    : bmds(mer_op_t::nb_mers)
    , nbmds(mer_op_t::nb_mers)
    , done(mer_op_t::nb_fmoves)
    , fmoves(mer_op_t::nb_fmoves)
    , nfmoves(mer_op_t::nb_fmoves)
    {}

    static bool has_fm(const std::vector<tristate_t> mds, mer_t fm) {
        for(mer_type b = 0; b < mer_op_t::alpha; ++b) {
            if(mds[mer_op_t::lc(fm, b)] != yes)
                return false;
        }
        return true;
    }

    // Given an MDS, fill bmds as a 1-hot vector representing mds and fill
    // fmoves as the equivalent order list of F-moves
    void mds2fmoves(const std::vector<mer_t>& mds) {
        std::fill(bmds.begin(), bmds.end(), nil);
        std::fill(done.begin(), done.end(), false);
        std::vector<mer_type> fms;
        for(auto m : mds) {
            bmds[m] = yes;
            const auto fm = mer_op_t::fmove(m);
            if(has_fm(bmds, fm))
                fms.push_back(fm);
        }

        mer_type i = 0;
        for( ; !fms.empty(); ++i) {

            mer_type fm = fms.back();
            fms.pop_back();
            fmoves[i] = fm;

            assert2(!done[fm], "Already done F-move " << fm);
            done[fm] = true;
            do_fmove(fm, bmds);

            const auto nm = mer_op_t::nmer(fm);
            for(mer_type b = 0; b < mer_op_t::alpha; ++b) {
                const auto nfm = mer_op_t::fmove(mer_op_t::rc(nm, b));
                if(!done[nfm] && has_fm(bmds, nfm))
                    fms.push_back(nfm);
            }
        }
        assert2(i == mer_op_t::nb_fmoves, "Too few F-moves done " << i);
    }


    static void do_fmove(mer_type fm, std::vector<tristate_t>& bmds) {
        for(mer_type b = 0; b < mer_op_t::alpha; ++b) {
            const auto m = mer_op_t::lc(fm, b);
            assert2(bmds[m] != no, "Left companion not present " << fm << ' ' << b);
            bmds[m] = no;
            bmds[mer_op_t::nmer(m)] = yes;
        }
    }

    // fromFmoves: initialize from a vector of fms. !!! Destructive !!!
    void fromFmoves(std::vector<mer_type>& fms) {
        assert2(fms.size() == mer_op_t::nb_fmoves, "Wrong number of f-moves given " << fms.size());
        fmoves.swap(fms);
    }

    // fromFmoves:: initialize from a vector of fms. !!! Destructive !!! (well,
    // this version is not, but assume fromFmoves is any)
    void fromFmoves(const std::vector<mer_t>& fms) {
        assert2(fms.size() == mer_op_t::nb_fmoves, "Wrong number of f-moves given " << fms.size());
        std::copy(fms.begin(), fms.end(), fmoves.begin());
    }

    // Given that fmoves is properly filled (e.g., by a call to
    // mds2fmoves() or fromFmoves(), and given a valid I-move imove (e.g., as found by
    // imoves_type.imoves()), create a 1-hot representation in nbmds of an MDS
    // in the component after traversing the I-move. nfmoves is a valid ordered
    // list of FMs in that new component as well.
    void traverse_imove(const imove_t& imove) {
        std::cout << "traverse " << imove << " | " << fmoves << std::endl;

        std::fill(bmds.begin(), bmds.end(), nil);
        std::fill(nbmds.begin(), nbmds.end(), nil);
        fm_listA.clear();
        nfmoves.clear();
        std::set<mer_type> targets;

        mask_t mask = 1;
        for(mer_t b = 0; b < mer_op_t::alpha; ++b, mask <<= 1) {
            if((imove.im & mask) != 0) {
                nbmds[mer_op_t::nmer(imove.fm, b)] = yes; // Mer known in new component
            } else {
                targets.insert(mer_op_t::lc(imove.fm, b));
            }
        }
        std::cout << "targets " << targets << std::endl;
        std::cout << "start nbmds " << nbmds << std::endl;

        const mer_type fmi = std::find(fmoves.cbegin(), fmoves.cend(), imove.fm) - fmoves.cbegin();
        nfmoves.push_back(imove.fm);
        mer_type i = 0;
        // Apply F-moves until all the targets (the nodes going around their
        // PCRs) are reached.
        for( ; i < mer_op_t::nb_fmoves && !targets.empty(); ++i) {
            mer_type j = fmi + i;
            if(j >= mer_op_t::nb_fmoves) j -= mer_op_t::nb_fmoves;
            const mer_type nfm = fmoves[j];
            std::cout << "nfm " << nfm << '\n';

            bool touch_nbmds = false;
            for(mer_type b = 0; !touch_nbmds && b < mer_op_t::alpha; ++b)
                touch_nbmds = nbmds[mer_op_t::lc(nfm, b)] == yes;
            if(touch_nbmds) {
                nfmoves.push_back(nfm);
                do_fmove(nfm, nbmds);
                std::cout << "\tListB\n";
            } else {
                fm_listA.push_back(nfm);
                do_fmove(nfm, bmds);
                for(auto it = targets.cbegin(); it != targets.cend(); ) {
                    if(bmds[*it] == yes) {
                        auto cit = it; // Have to copy first as .erase will invalidate the iterator
                        ++it;
                        targets.erase(cit);
                    } else {
                        ++it;
                    }
                }
                std::cout << "\tListA\n";
            }
        }

        // Continue doing the F-moves in order and apply to nbmds, then apply
        // from listA (except fm)
        for( ; i < mer_op_t::nb_fmoves; ++i) {
            mer_type j = fmi + i;
            if(j >= mer_op_t::nb_fmoves) j -= mer_op_t::nb_fmoves;
            const mer_type nfm = fmoves[j];
            do_fmove(nfm, nbmds);
            assert2(nfmoves.size() < mer_op_t::nb_fmoves, "Too many F-moves in new component, fmoves");
            nfmoves.push_back(nfm);
        }

        for(mer_type j = 1; j < fm_listA.size(); ++j) {
            const mer_type nfm = fm_listA[j];
            do_fmove(nfm, nbmds);
            assert2(nfmoves.size() < mer_op_t::nb_fmoves, "Too many F-moves in new component, listA");
            nfmoves.push_back(nfm);
        }

        // fm should be doable in nbmds
#ifndef NDEBUG
        std::cout << "nbmds " << nbmds << std::endl;
        for(mer_type b = 0; b < mer_op_t::alpha; ++b)
            assert2(nbmds[mer_op_t::lc(imove.fm, b)] == yes, "fm should be doable in new component " << imove.fm);
#endif
    }
};


#endif // MDS_OP_H_