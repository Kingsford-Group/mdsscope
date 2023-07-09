#ifndef MDS_OP_H_
#define MDS_OP_H_

#include "dbg.hpp"
#include "mer_op.hpp"
#include "common.hpp"

template<typename mer_op_type>
struct mds_op_type {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;

    std::vector<tristate_t> bmds;
    std::vector<bool> done;
    std::vector<mer_t> fmoves;

    mds_op_type()
    : bmds(mer_op_t::nb_mers)
    , done(mer_op_t::nb_fmoves)
    , fmoves(mer_op_t::nb_fmoves)
    {}

    static bool has_fm(const std::vector<tristate_t> mds, mer_t fm) {
        for(mer_type b = 0; has_fm && b < mer_op_t::alpha; ++b) {
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
};


#endif // MDS_OP_H_