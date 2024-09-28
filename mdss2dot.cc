#include <cstdlib>
#include <fstream>
#include <string>
#include <set>

#include "argparse.hpp"
#include "mer_op.hpp"
#include "common.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

struct MDSs2DOTArgs : argparse::Args {
    std::string& mds_arg = arg("DS list file");

    void welcome() override {
        std::cout << "Generate dot file from MDS list" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    const auto args = argparse::parse<MDSs2DOTArgs>(argc, argv);

    std::ifstream is(args.mds_arg);
    if(!is.good()) {
        std::cerr << "Failed to open " << args.mds_arg << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "digraph {\n";

    std::string line, token;
    std::set<mer_t> mds, dest;
    while(std::getline(is, line)) {
        mds.clear();

        // Read a line <-> MDS
        size_t start = 0, end;
        while((end = line.find(',', start)) != std::string::npos) {
            mds.insert(atoi(line.c_str() + start));
            start = end + 1;
        }
        mds.insert(atoi(line.c_str() + start));

        // Find potential F-move in MDS
        for(auto m : mds) {
            if(mer_ops::lb(m) != 0) continue;
            bool has_fmove = true;
            for(mer_t b = 1; has_fmove && b < mer_ops::alpha; ++b)
                has_fmove = mds.find(mer_ops::lc(m, b)) != mds.end();
            if(has_fmove) {
                dest = mds;
                for(mer_t b = 0; b < mer_ops::alpha; ++b) {
                    dest.erase(mer_ops::lc(m, b));
                    dest.insert(mer_ops::nmer(m, b));
                }
                std::cout << "  n" << join(mds, '_') << " -> n" << join(dest, '_') << ";\n";
            }
        }
    }

    std::cout << "}\n";

    return EXIT_SUCCESS;
}
