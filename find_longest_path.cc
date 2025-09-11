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
typedef mds_op_type<mer_ops> mds_ops;
typedef longest_path_type<mer_ops> longest_path;

struct LongestPathArgs : argparse::Args {
    bool& pairs = flag("p,pairs", "Print pairs of nodes of maximum path");
    bool& all_pair = flag("a,all", "Print the path start and end for all rfm");
    std::optional<const char*>& mds_arg = kwarg("f,mds", "File with MDS");
    std::vector<const char*>& comp_arg = arg("component").set_default("");

    void welcome() override {
        std::cout << "Find longest remaining path" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    const auto args = argparse::parse<LongestPathArgs>(argc, argv);
    std::vector<tristate_t> bmds;
	std::vector<mer_t> fms, rfms;
    longest_path lp;

    const auto mds = args.mds_arg ? mds_from_file<mer_t>(*args.mds_arg) : mds_from_arg<mer_t>(args.comp_arg);
    mds_ops::from_mds_fms(mds, bmds, fms, &rfms);
    const auto longest = lp.longest_path(bmds, fms);
    std::cout << (size_t)longest << '\n';

    if(args.all_pair || args.pairs) {
        for(const auto rfm : rfms) {
            for(mer_t b = 0; b < mer_ops::alpha; ++b) {
                const auto m = mer_ops::lc(rfm, b);
                if(args.all_pair || lp.paths[m] == longest)
                    std::cout << (size_t)lp.paths[m] << ": " << joinT<size_t>(*lp.origin[m], ',') << " -> " << (size_t)m << '\n';
            }
        }
    }

    return EXIT_SUCCESS;
}
