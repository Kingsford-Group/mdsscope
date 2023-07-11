#ifndef MDS_OP_H_
#define MDS_OP_H_

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
        std::fill(bmds.begin(), bmds.end(), nil);
        std::fill(nbmds.begin(), nbmds.end(), nil);
        fm_listA.clear();
        nfmoves.clear();

        mask_t mask = 1;
        for(mer_t b = 0; b < mer_op_t::alpha; ++b, mask <<= 1) {
            if(imove.im & mask == 0)
                nbmds[mer_op_t::nmer(imove.fm, b)] = yes; // Mer known in new component
        }


    }
};


#endif // MDS_OP_H_