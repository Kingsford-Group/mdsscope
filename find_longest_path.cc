#include "argparse.hpp"
#include <iostream>

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
#include "misc.hpp"
#include "longest_path.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef longest_path_type<mer_ops> longest_path;

struct LongestPathArgs : argparse::Args {
    std::optional<const char*>& mds_arg = kwarg("f,mds", "File with MDS");
    std::vector<const char*>& comp_arg = arg("component");

    void welcome() override {
        std::cout << "Find longest remaining path" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    const auto args = argparse::parse<LongestPathArgs>(argc, argv);
    longest_path lp;

    const auto mds = args.mds_arg ? mds_from_file<mer_t>(*args.mds_arg) : mds_from_arg<mer_t>(args.comp_arg);
    std::cout << (size_t)lp.longest_path(mds) << '\n';

    return EXIT_SUCCESS;
}
