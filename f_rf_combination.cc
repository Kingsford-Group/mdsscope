#include "argparse.hpp"
#include "common.hpp"
#include "connected_components.hpp"
#include "longest_path.hpp"
#include "mds_op.hpp"
#include "misc.hpp"
#include <cstdlib>

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef mds_op_type<mer_ops> mds_ops;
typedef connected_components_type<mer_ops> connected_components;
typedef longest_path_type<mer_ops> longest_path;

struct F_RF_CombintionArgs : argparse::Args {
   	std::optional<const char*>& mds_arg = kwarg("f,mds", "File with MDS");
    std::vector<const char*>& comp_arg = arg("component").set_default("");

    bool& debug = flag("d,debug", "Debug: print all combination found");
    bool& print = flag("p,print", "Print the combination if all fm in a combination");

    void welcome() override {
        std::cout << "Find F/RF combinations in the MDS" << std::endl;
    }
};

int main(int argc, char* argv[]) {
    const auto args = argparse::parse<F_RF_CombintionArgs>(argc, argv);
    std::vector<tristate_t> bmds;
	std::vector<mer_t> fms, rfms;
	connected_components cc;
	longest_path lp;

    const auto mds = args.mds_arg ? mds_from_file<mer_t>(*args.mds_arg) : mds_from_arg<mer_t>(args.comp_arg);
	mds_ops::from_mds_fms(mds, bmds, fms, &rfms);

	std::cout << "fms: " << joinT<size_t>(fms, ',') << '\n'
              << "rfms: " << joinT<size_t>(rfms, ',') << '\n';

	std::set<mer_t> rfmss(rfms.begin(), rfms.end());
	// For each lc of fm, check if the node is also part of a rfm.
	bool all_fm_rfms = true; // All fm part of a combination
	for(const auto fm : fms) {
	    bool this_fm_rfms = false;
		mer_t rfm;
	    for(mer_t b = 0; b < mer_ops::alpha; ++b) {
			rfm = mer_ops::rfmove(mer_ops::lc(fm, b));
			if(rfmss.contains(rfm)) {
			    this_fm_rfms = true;
				break;
			}
		}
	    if(args.debug) {
			if(this_fm_rfms)
    		    std::cout << "Combo: " << (size_t)fm << ' ' << (size_t)rfm << ' ' << all_fm_rfms << '\n';
			else
			    std::cout << "Not: " << (size_t)fm << '\n';
		}
		all_fm_rfms = all_fm_rfms && this_fm_rfms;
	}

	if(all_fm_rfms) {
	  std::cout << "All fm in combo" << std::endl;
	}

	cc.compute_components(bmds);
	// Not very efficient way to display the components, but whatever
	std::set<mer_t> comps;
	// Collect the component indices
	for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
		if(bmds[m] == yes) continue;
		comps.insert(cc.components.find(m));
	}
	std::cout << "Components " << comps.size() << ":\n";

	// Output the components one by one
	for(const auto idx : comps) {
		std::cout << (size_t)idx << ':';
		for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
			if(cc.components.find(m) == idx && bmds[m] != yes)
				std::cout << ' ' << (size_t)m;
		}
		std::cout << '\n';
	}

	std::cout << "Max length " << (size_t)lp.longest_path(bmds, fms) << '\n';
	for(const auto rfm : rfms) {
        for(mer_t b = 0; b < mer_ops::alpha; ++b) {
            const auto m = mer_ops::lc(rfm, b);
            if(bmds[m] == yes) { // This should only happen for the homopolymer RFM
                assert2(mer_ops::is_homopolymer_fm(rfm), "Left companion in set for non homopolymer RFM");
                continue;
            }
            assert2(lp.origin[m] != lp.origin_set.end(), "Origin not set");
            std::cout << (size_t)lp.paths[m] << ':';
           for(const auto start : *lp.origin[m])
               std::cout << ' ' << (size_t)start << '(' << (size_t)mer_ops::rfmove(start) << ')';
           std::cout << " -> " << (size_t)m << '(' << (size_t)mer_ops::fmove(m) << ")\n";
        }
    }

	return EXIT_SUCCESS;
}
