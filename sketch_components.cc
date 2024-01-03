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
	mer_t current, scc_index;
	std::vector<mer_t> index;
	std::vector<mer_t> lowlink;
	std::vector<bool> onstack;

	tarjan_scc()
	: index(mer_ops::nb_mers, undefined)
	, lowlink(mer_ops::nb_mers, undefined)
	, onstack(mer_ops::nb_mers, false)
		{}

	// Find strongly connected component in the de Bruijn graph minus a set. The
	// set is given by the indicator function `fn`. For each component, call
	// `new_scc(scc_index)` and then call `new_node(m)` for all the mers in the
	// component.
	template<typename Fn, typename E1, typename E2>
	void scc_iterate(Fn fn, E1 new_scc, E2 new_node) {
		std::fill(index.begin(), index.end(), undefined);
		std::fill(lowlink.begin(), lowlink.end(), undefined);
		std::fill(onstack.begin(), onstack.end(), false);

		current = 0;
		scc_index= 0;
		stack.clear();

		for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
			if(index[m] == undefined && !fn(m))
				strong_connect(fn, m, new_scc, new_node);
		}
	}

	template<typename Fn>
	std::vector<std::vector<mer_t>> scc_components(Fn fn) {
		std::vector<std::vector<mer_t>> res;
		auto new_scc = [&res](mer_t index) { res.emplace_back(); };
		auto new_node = [&res](mer_t m) { res.back().push_back(m); };
		scc_iterate(fn, new_scc, new_node);
		return res;
	}

	template<typename Fn>
	std::pair<mer_t, mer_t> scc_counts(Fn fn) {
		mer_t nb_scc = 0, nb_mers = 0;
		auto new_scc = [&nb_scc](mer_t index) { ++nb_scc; };
		auto new_mer = [&nb_mers](mer_t m) { ++nb_mers; };
		scc_iterate(fn, new_scc, new_mer);
		return std::make_pair(nb_scc, nb_mers);
	}

	template<typename Fn>
	std::pair<mer_t, mer_t> scc_append(Fn fn, std::vector<mer_t>& mers) {
		mer_t nb_scc = 0, nb_mers = 0;
		auto new_scc = [&nb_scc](mer_t index) { ++nb_scc; };
		auto new_mer = [&](mer_t m) { ++nb_mers; mers.push_back(m); };
		scc_iterate(fn, new_scc, new_mer);
		return std::make_pair(nb_scc, nb_mers);
	}

private:
	template<typename Fn, typename E1, typename E2>
	void strong_connect(Fn fn, mer_t m, E1 new_scc, E2 new_node) {
		index[m] = current;
		lowlink[m] = current;
		++current;
		stack.push_back(m);
		onstack[m] = true;

		for(mer_t b = 0; b < mer_ops::alpha; ++b) {
			mer_t nmer = mer_ops::nmer(m, b);
			if(fn(nmer)) continue;

			if(index[nmer] == undefined) {
				strong_connect(fn, nmer, new_scc, new_node);
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

struct mer_or_rc_in_set {
	const std::unordered_set<mer_t>& set;
	mer_or_rc_in_set(const std::unordered_set<mer_t>& s) : set(s) {}
	bool operator()(mer_t m) const { return set.find(m) != set.cend() || set.find(mer_ops::reverse_comp(m)) != set.cend(); }
};

// // Transform a set into a "canonical set" based on whether PCRs are contained
// // within the k-nonical space or not.
// std::unordered_set<mer_t> mixed_pcr_transform(const std::unordered_set<mer_t>& set) {
// 	std::cout << "input " << joinT<size_t>(set, ',') << '\n';
// 	std::unordered_set<mer_t> res;

// 	for(mer_t m : set) {
// 		mer_t rcm = mer_ops::reverse_comp(m);
// 		unsigned char pcr_type = 0; // first bit: in k-nonical, second bit: in rc
// 		pcr_type |= (m < rcm) ? 0x1 : ((m > rcm) ? 0x2 : 0x3) ;

// 		mer_t cur = mer_ops::nmer(m);
// 		mer_t rcc = mer_ops::pmer(rcm);
// 		for( ; cur != m && pcr_type != 0x3; cur = mer_ops::nmer(cur), rcc = mer_ops::pmer(rcc)) {
// 			pcr_type |= (cur < rcc) ? 0x1 : ((cur > rcc) ? 0x2 : 0x3);
// 			std::cout << '\t' << (size_t)cur << ' ' << (size_t)rcc << ' ' << (size_t)pcr_type << '\n';
// 		}
// 		std::cout << "mer " << (size_t)m << ' ' << (size_t)rcm << ' ' << (size_t)pcr_type << '\n';
// 		switch(pcr_type) {
// 		case 0x1: // All in k-nonical
// 			std::cout << "insert " << (size_t)m << '\n';
// 			res.insert(m);
// 			break; //
// 		case 0x3: // Mixed PCR
// 			std::cout << "insert " << (size_t)std::min(m, rcm) << '\n';
// 			res.insert(std::min(m, rcm));
// 			break;
// 		default:
// 			// All in rc: do nothing
// 			std::cout << "nothing\n";
// 			break;
// 		}
// 	}
// 	std::cout << "res " << joinT<size_t>(res, ',') << '\n';
// 	return res;
// }

// Symmetrize a set. Well defined from MDSs. Not so sure about other sets
std::unordered_set<mer_t> symmetrize(const std::unordered_set<mer_t>& set) {
	std::unordered_set<mer_t> res;



	return res;
}

int main(int argc, char* argv[]) {
	sketch_components args(argc, argv);
	const auto mer_set = read_set(args);

	const is_in_set set_fn(mer_set.first);
	const can_is_in_set can_fn(mer_set.first);
	const mer_or_rc_in_set or_fn(mer_set.first);

	tarjan_scc comp_scc;
	const auto res = comp_scc.scc_counts(set_fn);
	const auto can_res = comp_scc.scc_counts(can_fn);
	const auto or_res = comp_scc.scc_counts(or_fn);

	std::cout << "f\t" << (size_t)res.first << ',' << (size_t)res.second << '\n'
			  << "fc\t" << (size_t)can_res.first << ',' << (size_t)can_res.second << '\n'
			  << "f + frc\t" << (size_t)or_res.first << ',' << (size_t)or_res.second << '\n';

	return EXIT_SUCCESS;
}
