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

typedef mer_op_type<K, ALPHA> mer_ops;


int main(int argc, char* argv[]) {
    show_backtrace();
    std::ios::sync_with_stdio(false);
    traverse_comp args(argc, argv);

    imoves_type<mer_ops> imoves_op;
    mds_op_type<mer_ops> mds_op;

    queue_elt<mer_ops> current;
    comp_queue<mer_ops> queue(args.comps_arg);

    const auto start(mds_from_arg<mer_type>(args.comp_arg));
    current.index = 0;
    current.ims = imoves_op.imoves(start);
    mds_op.mds2fmoves(start);
    current.fms = mds_op.fmoves;

    queue.enqueue(current);

    while(!queue.empty()) {
        if(!queue.dequeue(current))
            throw std::runtime_error("Failed to dequeue element");
        std::cout << current << '\n';
        mds_op.fromFmoves(current.fms);

    }

    return 0;
}