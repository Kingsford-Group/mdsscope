// For a particular sketching method than may not define an MDS, output the
// strongly connected components of the de Bruijn graph minus the set of
// selected k-mers.

#include <unordered_set>
#include <vector>
#include <stack>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cstddef> // for std::ptrdiff_t
#include "misc.hpp"
#include "common.hpp"
#include "sketch_components.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;


std::pair<std::unordered_set<mer_t>, mer_t> read_set(const sketch_components& args) {
	const auto mers = args.set_given ? read_mds_from_file<mer_t>(args.set_arg) : mds_from_arg<mer_t>(args.comp_arg);
	std::unordered_set<mer_t> set;
	mer_t size = 0;
	for(const auto m : mers) {
		if(set.insert(m).second) ++size;
	}
	return std::make_pair(set, size);
}

// Function checking if a mer is in the set of selecting mers.
typedef bool (*in_set_fn)(mer_t);

struct tarjan_scc {
	static const mer_t undefined = mer_ops::nb_mers;
	std::vector<mer_t> stack;
	mer_t current;
	std::vector<mer_t> index;
	std::vector<mer_t> lowlink;
	std::vector<bool> onstack;

	tarjan_scc()
	: index(mer_ops::nb_mers, undefined)
	, lowlink(mer_ops::nb_mers, undefined)
	, onstack(mer_ops::nb_mers, false)
		{}

	template<typename Fn>
	std::vector<std::vector<mer_t>> scc(Fn fn) {
		std::vector<std::vector<mer_t>> res;
		std::fill(index.begin(), index.end(), undefined);
		std::fill(lowlink.begin(), lowlink.end(), undefined);
		std::fill(onstack.begin(), onstack.end(), false);

		current = 0;
		stack.clear();

		for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
			if(index[m] == undefined && !fn(m))
				strong_connect(fn, m, res);
		}
		return res;
	}

	template<typename Fn>
	void strong_connect(Fn fn, mer_t m, std::vector<std::vector<mer_t>>& res) {
		index[m] = current;
		lowlink[m] = current;
		++current;
		stack.push_back(m);
		onstack[m] = true;

		for(mer_t b = 0; b < mer_ops::alpha; ++b) {
			mer_t nmer = mer_ops::nmer(m, b);
			if(fn(nmer)) continue;

			if(index[nmer] == undefined) {
				strong_connect(fn, nmer, res);
				lowlink[m] = std::min(lowlink[m], lowlink[nmer]);
			} else if(onstack[nmer]) {
				lowlink[m] = std::min(lowlink[m], index[nmer]);
			}
		}

		if(lowlink[m] == index[m]) {
			// A single node is not an SCC, unless has a self loop (homopolymers)
			if(stack.back() == m && !mer_ops::is_homopolymer(m)) {
				onstack[m] = false;
				stack.pop_back();
			} else {
				std::vector<mer_t> new_scc;
				while(true) {
					const auto mm = stack.back();
					stack.pop_back();
					onstack[mm] = false;
					new_scc.push_back(mm);
					if(mm == m) break;
				}
				res.push_back(std::move(new_scc));
			}
		}
	}
};

struct is_in_set {
	const std::unordered_set<mer_t>& set;
	is_in_set(const std::unordered_set<mer_t>& s) : set(s) {}
	bool operator()(mer_t m) const { return set.find(m) != set.cend(); }
};

int main(int argc, char* argv[]) {
	sketch_components args(argc, argv);
	const auto mer_set = read_set(args);

	const is_in_set set_fn(mer_set.first);

	tarjan_scc comp_scc;
	const auto res = comp_scc.scc(set_fn);

	mer_t index = 0;
	for(const auto& scc : res) {
		std::cout << (size_t)(index++) << ' ' << scc.size() << ": " << joinT<size_t>(scc, ',') << '\n';
	}

	return EXIT_SUCCESS;
}
