// For a particular sketching method than may not define an MDS, output the
// strongly connected components of the de Bruijn graph minus the set of
// selected k-mers.

#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <vector>
#include <stack>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>

#include "misc.hpp"
#include "common.hpp"
#include "sketch_components.hpp"
#include "tarjan_scc.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

// Function checking if a mer is in the set of selecting mers.
typedef bool (*in_set_fn)(mer_t);

struct is_in_set {
	const std::unordered_set<mer_t>& set;
	is_in_set(const std::unordered_set<mer_t>& s) : set(s) {}
	bool operator()(mer_t m) const { return set.find(m) != set.cend(); }
};

struct can_is_in_set {
	const std::unordered_set<mer_t>& set;
	can_is_in_set(const std::unordered_set<mer_t>& s) : set(s) {}
	bool operator()(mer_t m) const { return set.find(mer_ops::canonical(m)) != set.cend(); }
};

struct is_in_union {
	const std::unordered_set<mer_t>& set;
	is_in_union(const std::unordered_set<mer_t>& s) : set(s) {}
	bool operator()(mer_t m) const { return set.find(m) != set.cend() || set.find(mer_ops::reverse_comp(m)) != set.cend(); }
};

int main(int argc, char* argv[]) {
	sketch_components args(argc, argv);
	const auto mer_set = get_mds<std::unordered_set<mer_t>>(args.sketch_file_arg, args.sketch_arg);

	mer_t components = 0, in_components = 0, visited = 0;
	size_t updates = 0;
	auto progress = [&]() {
		if(args.progress_flag) {
			if(updates % 1024 == 0) {
				std::cerr << '\r'
						  << "comps " << std::setw(10) << components
						  << " in_comps " << std::setw(10) << in_components
						  << " visited " << std::setw(10) << visited
						  << ' ' << std::setw(5) << std::fixed << std::setprecision(1) << (100.0 * (double)visited / mer_ops::nb_mers) << '%'
						  << std::flush;
			}
			++updates;
		}
	};
	auto new_scc = [&components,&progress](mer_t m) { ++components; progress(); };
	auto new_node = [&in_components,&progress](mer_t m) { ++in_components; progress(); };
	auto new_visit = [&visited,&progress](mer_t m) { ++visited; progress(); };

	tarjan_scc<mer_ops> comp_scc;
	if(args.canonical_flag) {
		const can_is_in_set can_fn(mer_set);
		comp_scc.scc_iterate(can_fn, new_scc, new_node, new_visit);
	} else if(args.union_flag) {
		const is_in_union union_fn(mer_set);
		comp_scc.scc_iterate(union_fn, new_scc, new_node, new_visit);
	} else {
		const is_in_set set_fn(mer_set);
		comp_scc.scc_iterate(set_fn, new_scc, new_node, new_visit);
	}

	if(args.progress_flag)
		std::cerr << std::endl;

	std::cout << (size_t)components << ',' << (size_t)in_components << '\n';

	return EXIT_SUCCESS;
}
