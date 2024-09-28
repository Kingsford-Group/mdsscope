#include <iostream>
#include <string>
#include <vector>

#include "argparse.hpp"
#include "mer_op.hpp"
#include "mds_op.hpp"
#include "misc.hpp"
#include "common.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

struct FMS2MDSArgs : argparse::Args {

    void welcome() {
        std::cout << "Transform list of FMs to an mers in an MDS" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    const auto args = argparse::parse<FMS2MDSArgs>(argc, argv);
    mds_op_type<mer_ops> mds_op;

    std::string line;
    std::vector<const char*> l(1, nullptr);
    while(std::getline(std::cin, line)) {
        l[0] = line.c_str();
        const auto fms = mds_from_arg<mer_t>(l, false);
        mds_op.fmoves2mds(fms);
        std::cout << mds_op.bmds << '\n';
    }

    return 0;
}
