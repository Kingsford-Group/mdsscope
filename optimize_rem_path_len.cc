#include <cstdlib>
#include <random>

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
#include "mds_op.hpp"
#include "imoves.hpp"
#include "imove_signature.hpp"
#include "longest_path.hpp"
#include "file_queue.hpp"
#include "misc.hpp"

#include "optimize_rem_path_len.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

template<typename mer_op_type>
struct element {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;
    typedef imove_sig_type<mer_op_type> imove_sig_t;
    typedef mds_op_type<mer_op_type> mds_t;

    mer_t path_len;
    std::vector<tristate_t> bmds;
    std::vector<mer_t> fms; // All f-moves equivalent to bmds
    std::vector<mer_t> fmoves; // possible f-moves
    imove_sig_t ims;
};

volatile bool interrupt = false;
void int_handler(int sig) {
    interrupt = true;
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    imoves_type<mer_ops> imoves_op;
    mds_op_type<mer_ops> mds_op;
    longest_path_type<mer_ops> longest_path;

    std::mt19937_64 rand_gen((std::random_device())());
    std::uniform_real_distribution<float> rand_unit(0.0, 1.0);

    optimize_rem_path_len args(argc, argv);

    element<mer_ops> current, best;
    std::vector<mer_t> fms;

    const auto start(mds_from_arg<mer_type>(args.comp_arg));
    mds_op.from_mds_fms(start, current.bmds, current.fms);
    mds_op.mds2fmoves(start);
    current.fmoves = mds_op.fmoves;
    current.path_len = longest_path.longest_path(current.bmds, current.fms);
    current.ims = imoves_op.imoves(start);
    best = current;

    std::cout << current.path_len << std::endl;

    // Do simulated annealing. Energy is length of remaining path. T starts at 1
    // and updated by T_{n+1} = \lambda*T_n. DE = E_{new} - E_{current}.
    // Probability of acceptence is 1 if DE < 0 and exp(-DE/T))/2 if DE >= 0.
    double temp = 1.0;

    // Do one step of I-move and lower temperature. Within each of this steps,
    // do at most 2*K F-moves without lowering the temperature.
    for(uint64_t iteration = 0; iteration < args.iteration_arg; ++iteration, temp *= args.lambda_arg) {
        std::uniform_int_distribution<int> rand_im(0, current.ims.size() - 1);
        const auto im = current.ims[rand_im(rand_gen)];
        mds_op.fromFmoves((const std::vector<mer_t>&)current.fmoves);
        mds_op.traverse_imove(im);
        mds_op.from_bmds_fms(mds_op.nbmds, fms);

        mer_t nlp = longest_path.longest_path(mds_op.nbmds, fms);
        // std::cout << "im " << im << ' ' << nlp << " | " << current.bmds << " | " << mds_op.nbmds << std::endl;
        // Do a some numbers of F-move to find a low longest path within the
        // component
        for(unsigned int fmi = 0; fmi < 2 * mer_ops::k; ++fmi) {
            const mer_t fm = mds_op.nfmoves[fmi];
            mds_op.do_fmove(fm, mds_op.nbmds);
            fms.erase(std::find(fms.begin(), fms.end(), fm));
            for(mer_t b = 0; b < mer_ops::alpha; ++b) {
                mer_t nfm = mer_ops::fmove(mer_ops::nmer(fm, b));
                if(mds_op.has_fm(mds_op.nbmds, nfm))
                   fms.push_back(nfm);
            }
            mer_t nnlp = longest_path.longest_path(mds_op.nbmds, fms);
            nlp = std::min(nlp, nnlp);
        }

        if(nlp >= current.path_len && rand_unit(rand_gen) > std::exp(- (float)(nlp - current.path_len) / temp))
            continue;

        // Update current
        current.path_len = nlp;
        current.bmds = mds_op.nbmds;
        current.fmoves = mds_op.nfmoves;
        current.fms = fms;
        current.ims = imoves_op.imoves(current.bmds);

        if(current.path_len < best.path_len)
            best = current;
    }

    std::cout << best.path_len << std::endl;
    return EXIT_SUCCESS;
}
