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
#include "misc.hpp"
#include "backtrace.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;


int main(int argc, char* argv[]) {
    show_backtrace();
    std::ios::sync_with_stdio(false);
    traverse_comp args(argc, argv);

    const auto start(mds_from_arg<mer_type>(args.comp_arg));
    imoves_type<mer_ops> imoves_op;

    const auto ims = imoves_op.imoves(start);
    std::cout << start << std::endl;
    std::cout << ims << std::endl;


    return 0;
}