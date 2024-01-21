#include <iostream>
#include <unordered_set>

#include "misc.hpp"
#include "common.hpp"
#include "tarjan_scc.hpp"

#include "opt_canon.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;


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

struct is_in_opt {
	const std::unordered_set<mer_t>& set1;
	const std::unordered_set<mer_t>& set2;
	is_in_opt(const std::unordered_set<mer_t>& s1, const std::unordered_set<mer_t>& s2)
		: set1(s1)
		, set2(s2)
		{}
	bool operator()(mer_t m) const {
		return (set1.find(m) != set1.cend() || set1.find(mer_ops::reverse_comp(m)) != set1.cend()) && \
			(set2.find(m) == set2.end());
	}
};

int main(int argc, char* argv[]) {
	opt_canon args(argc, argv);
	const auto mer_set = get_mds<std::unordered_set<mer_t>>(args.sketch_file_arg, args.sketch_arg);

	tarjan_scc<mer_ops> comp_scc;
	const is_in_set set_fn(mer_set);
	const can_is_in_set can_fn(mer_set);
	const auto counts_orig = comp_scc.scc_counts(set_fn);
	std::cout << counts_orig.first << ' ' << counts_orig.second << '\n';

	std::unordered_set<mer_t> remove;

	// Greedily create the remove set. Should compare the SCCs themselves, not
	// just their number and size. Sufficient for now....
	for(const auto m : mer_set) {
		if(can_fn(m)) continue; // Must be a super-set of canonical
		// Try adding m to remove. If increase SCCs, don't keep it
		remove.insert(m);
		const is_in_opt opt(mer_set, remove);
		const auto counts = comp_scc.scc_counts(opt);
		if(counts.first > counts_orig.first || counts.second > counts_orig.second)
			remove.erase(m);
	}

	std::cout << "orig   " << mer_set.size() << ": " << joinT<size_t>(mer_set, ',') << '\n';
	std::cout << "remove " << remove.size() << ": " << joinT<size_t>(remove, ',') << '\n';

	return EXIT_SUCCESS;
}
