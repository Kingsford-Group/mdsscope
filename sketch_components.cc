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
#include "progressbar.hpp"

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

struct tarjan_scc {
	static const mer_t undefined = mer_ops::nb_mers;
	std::vector<mer_t> stack;
	mer_t current, scc_index;
	std::vector<mer_t> index;
	std::vector<mer_t> lowlink;
	std::vector<bool> onstack;
	std::vector<std::pair<mer_t,unsigned>> callstack; // To linearize algorithm

	tarjan_scc()
	: index(mer_ops::nb_mers, undefined)
	, lowlink(mer_ops::nb_mers, undefined)
	, onstack(mer_ops::nb_mers, false)
		{}

	static void noprogress(mer_t m) { }

	// Find strongly connected component in the de Bruijn graph minus a set. The
	// set is given by the indicator function `fn`. For each component, call
	// `new_scc(scc_index)` and then call `new_node(m)` for all the mers in the
	// component.
	template<typename Fn, typename E1, typename E2, typename E3>
	void scc_iterate(Fn fn, E1 new_scc, E2 new_node, E3 new_visit = noprogress) {
		std::fill(index.begin(), index.end(), undefined);
		std::fill(lowlink.begin(), lowlink.end(), undefined);
		std::fill(onstack.begin(), onstack.end(), false);

		current = 0;
		scc_index= 0;
		stack.clear();

		for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
			if(index[m] == undefined && !fn(m))
				strong_connect(fn, m, new_scc, new_node, new_visit);
		}
	}

	template<typename Fn, typename E3>
	std::vector<std::vector<mer_t>> scc_components(Fn fn, E3 new_visit = noprogress) {
		std::vector<std::vector<mer_t>> res;
		auto new_scc = [&res](mer_t index) { res.emplace_back(); };
		auto new_node = [&res](mer_t m) { res.back().push_back(m); };
		scc_iterate(fn, new_scc, new_node, new_visit);
		return res;
	}

	template<typename Fn, typename E3>
	std::pair<mer_t, mer_t> scc_counts(Fn fn, E3 new_visit = noprogress) {
		mer_t nb_scc = 0, nb_mers = 0;
		auto new_scc = [&nb_scc](mer_t index) { ++nb_scc; };
		auto new_mer = [&nb_mers](mer_t m) { ++nb_mers; };
		scc_iterate(fn, new_scc, new_mer, new_visit);
		return std::make_pair(nb_scc, nb_mers);
	}

	template<typename Fn, typename E3>
	std::pair<mer_t, mer_t> scc_append(Fn fn, std::vector<mer_t>& mers, E3 new_visit = noprogress) {
		mer_t nb_scc = 0, nb_mers = 0;
		auto new_scc = [&nb_scc](mer_t index) { ++nb_scc; };
		auto new_mer = [&](mer_t m) { ++nb_mers; mers.push_back(m); };
		scc_iterate(fn, new_scc, new_mer, new_visit);
		return std::make_pair(nb_scc, nb_mers);
	}

private:
	// Non-recursive implementation of Tarjan algorithm to find SCCs. Calls
	// new_scc upon finding a new SCC, and then calls new_node for each node in
	// that SCC.
	template<typename Fn, typename E1, typename E2, typename E3>
	inline void strong_connect(Fn fn, mer_t m, E1 new_scc, E2 new_node, E3 new_visit) {
		new_visit(m);
		index[m] = current;
		lowlink[m] = current;
		++current;
		stack.push_back(m);
		onstack[m] = true;
		callstack.emplace_back(m, 0);

		unsigned b;
		while(!callstack.empty()) {
			std::tie(m, b) = callstack.back();
			if(b < mer_ops::alpha) {
				// Still edges to explore
				++callstack.back().second;
				const mer_t nmer = mer_ops::nmer(m, b);
				if(fn(nmer)) continue;

				if(index[nmer] == undefined) {
					// Explore neighbor: push stack
					new_visit(m);
					index[nmer] = current;
					lowlink[nmer] = current;
					++current;
					stack.push_back(nmer);
					onstack[nmer] = true;
					callstack.emplace_back(nmer, 0);
				} else if(onstack[nmer]) {
					lowlink[m] = std::min(lowlink[m], index[nmer]);
				}
			} else {
				if(lowlink[m] == index[m]) {
					// A single node is not an SCC, unless has a self loop (homopolymers)
					if(stack.back() == m && !mer_ops::is_homopolymer(m)) {
						onstack[m] = false;
						stack.pop_back();
					} else {
						new_scc(scc_index++);
						while(true) {
							const auto mm = stack.back();
							stack.pop_back();
							onstack[mm] = false;
							new_node(mm);
							if(mm == m) break;
						}
					}
				}
				// Pop callstack
				callstack.pop_back();
				if(!callstack.empty()) {
					std::tie(m, b) = callstack.back();
					mer_t nmer = mer_ops::nmer(m, b-1); // Value of nmer when pushed to callstack
					lowlink[m] = std::min(lowlink[m], lowlink[nmer]);
				}
			}
		}
	}
};
const mer_t tarjan_scc::undefined;

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

	tarjan_scc comp_scc;
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

	std::cout << components << ',' << in_components << '\n';

	return EXIT_SUCCESS;
}
