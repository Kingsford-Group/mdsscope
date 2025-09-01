#include "argparse.hpp"
#include "common.hpp"
#include "mds_op.hpp"
#include "misc.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "connected_components.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef mds_op_type<mer_ops> mds_ops;
typedef connected_components_type<mer_ops> connected_components;

struct ConnectedComponentsArgs : argparse::Args {
	std::optional<const char*>& mds_arg = kwarg("f,mds", "File with MDS");
	bool& print_comps = flag("p,print", "Print the components");
    std::vector<const char*>& comp_arg = arg("component").set_default("");

    void welcome() override {
        std::cout << "Find longest remaining path" << std::endl;
    }
};

int main(int argc, char* argv[]) {
	const auto args = argparse::parse<ConnectedComponentsArgs>(argc, argv);
    connected_components cc;
	std::vector<tristate_t> bmds;
	std::vector<mer_t> fms;

    const auto mds = args.mds_arg ? mds_from_file<mer_t>(*args.mds_arg) : mds_from_arg<mer_t>(args.comp_arg);
	mds_ops::from_mds_fms(mds, bmds, fms);
	cc.compute_components(bmds);
    std::cout << ((size_t)cc.components.nb_sets - mds.size()) << ' '
			  << (size_t)cc.components.nb_sets << ' '
			  << mds.size() << '\n';

	if(args.print_comps) {
		// Not very efficient way to display the components, but whatever
		std::set<mer_t> comps;
		// Collect the component indices
		for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
			if(bmds[m] == yes) continue;
			comps.insert(cc.components.find(m));
		}

		// Output the components one by one
		for(const auto idx : comps) {
			std::cout << (size_t)idx << ':';
			for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
				if(cc.components.find(m) == idx && bmds[m] != yes)
					std::cout << ' ' << (size_t)m;
			}
			std::cout << '\n';
		}
	}

    return EXIT_SUCCESS;
}
