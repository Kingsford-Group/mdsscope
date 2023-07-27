#include "traverse_comp.hpp"

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
#include "file_queue.hpp"
#include "misc.hpp"
#include "backtrace.hpp"

// Only used for debugging.
#ifndef NDEBUG
#include "pcr_info.hpp"
#endif

typedef mer_op_type<K, ALPHA> mer_ops;


int main(int argc, char* argv[]) {
    show_backtrace();
    std::ios::sync_with_stdio(false);
    traverse_comp args(argc, argv);

    imoves_type<mer_ops> imoves_op;
    mds_op_type<mer_ops> mds_op;

    queue_elt<mer_ops> current, nelt;
    comp_queue<mer_ops> queue(args.comps_arg);
    signatures_type<mer_ops> signatures;

#ifndef NDEBUG
    pcr_info_type<mer_ops> pcr_info;
#endif

    const auto start(mds_from_arg<mer_type>(args.comp_arg));
    assert2(pcr_info.check_mds(start), "Invalid starting MDS");
    // std::cout << "start " << start << std::endl;
    current.index = 0;
    current.ims = imoves_op.imoves(start);
    signatures.insert(std::make_pair(current.ims, current.index));
    // std::cout << "signatures " << signatures << std::endl;

    // std::cout << current.ims.size() << ": " << current.ims << std::endl;
    mds_op.mds2fmoves(start);
    current.fms = mds_op.fmoves;

    queue.enqueue(current);
    nelt.fms.resize(mer_ops::nb_mers);

    while(!queue.empty()) {
        if(!queue.dequeue(current))
            throw std::runtime_error("Failed to dequeue element");
        // std::cout << "current " << current << std::endl;
        mds_op.fromFmoves(current.fms);

        for(const auto im : current.ims) {
            // std::cout << "i-move " << im << std::endl;
            mds_op.traverse_imove(im);
            assert2(pcr_info.check_bmds(mds_op.nbmds), "Invalid bmds after traversing I-move");
            imoves_op.imoves(mds_op.nbmds, nelt.ims);
            // std::cout << "nfmoves " << mds_op.nfmoves << " imoves " << nelt.ims << std::endl;
            nelt.fms.swap(mds_op.nfmoves);
            const auto ires = signatures.insert(std::make_pair(nelt.ims, signatures.size()));
            // std::cout << "insert " << ires.second << ' ' << ires.first->second << " | " << ires.first->first << std::endl;
            // std::cout << "signatures " << signatures << std::endl;
            nelt.index = ires.first->second;
            // TODO: output dot file
            if(!ires.second) continue;
            queue.enqueue(nelt);
        }
    }

    return 0;
}
