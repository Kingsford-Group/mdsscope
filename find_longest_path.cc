#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "find_longest_path.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
#include "mds_op.hpp"
#include "misc.hpp"
#include "longest_path.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef longest_path_type<mer_ops> longest_path;

int main(int argc, char* argv[]) {
    find_longest_path args(argc, argv);
    longest_path lp;

    const auto mds = args.mds_given ? read_mds_from_file<mer_t>(args.mds_arg) : mds_from_arg<mer_t>(args.comp_arg);
    std::cout << (size_t)lp.longest_path(mds) << '\n';

    return EXIT_SUCCESS;
}
