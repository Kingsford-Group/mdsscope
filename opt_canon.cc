#include <iostream>
#include <unordered_set>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>

#include "misc.hpp"
#include "common.hpp"
#include "tarjan_scc.hpp"
#include "random_seed.hpp"

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

// True if mer m or its reverse complement is in s1, but neither are in s2.
struct is_in_opt {
	const std::unordered_set<mer_t>& set1;
	const std::unordered_set<mer_t>& set2;
	is_in_opt(const std::unordered_set<mer_t>& s1, const std::unordered_set<mer_t>& s2)
		: set1(s1)
		, set2(s2)
		{}
	bool operator()(mer_t m) const {
		const auto rcm = mer_ops::reverse_comp(m);
		return (set1.find(m) != set1.cend() || set1.find(rcm) != set1.cend()) && \
			(set2.find(m) == set2.end() && set2.find(rcm) == set2.end());
	}
};

template<typename R, typename C>
R to_canonical_set(const C& mers) {
	R res;

	for(const auto m : mers) {
		const auto& rc = mer_ops::reverse_comp(m);
		if(rc < m) continue;
		res.insert(m);
		res.insert(rc);
	}

	return res;
}

template<typename C>
size_t canonicalize_size(const C& mers) {
	size_t size = 0;
	for(const auto m : mers) {
		const auto rcm = mer_ops::reverse_comp(m);
		// Add two for canonical k-mers (1 for itself, 1 for its rc), unless it is self rc (then add only 1)
		if(m <= rcm) {
			++size;
			if(m < rcm)
				++size;
		}
	}

	return size;
}

template<typename C, typename R>
size_t union_size(const C& mers, const R* remove) {
	size_t size = 0;
	for(const auto m : mers) {
		if(remove && remove->find(m) != remove->end()) continue;
		++size; // Add one for the k-mer itself
		const auto rcm = mer_ops::reverse_comp(m);
		if(m != rcm && mers.find(rcm) == mers.end())
			++size; // Add one for its rc if not in set
	}
	return size;
}

// Does a BFS to detect a new cycle in the de Bruijn graph minus a set. Starts
// from m and the reverse complement of m (rcm) and check for a loop back to m
// or back to rcm.
struct symm_bfs {
	std::queue<mer_t> queue;
	std::vector<bool> visited;

	static void noprogress(mer_t i) {}

	symm_bfs() : visited(mer_ops::nb_mers) {}

	template<typename Fn>
	bool has_cycle(Fn in_set, mer_t m) {
		std::cout << "Has cycle " << (size_t)m << std::endl;
		const auto rcm = mer_ops::reverse_comp(m);
		std::cout << "Clear" << std::endl;
		std::fill(visited.begin(), visited.end(), false);
		clear_queue();

		size_t nb_visited = 0;
		queue.push(m);
		while(true) { // Repeat with rcm eventually
			while(!queue.empty()) {
				if(nb_visited % 1000 == 0)
					std::cout << '\r' << nb_visited << std::flush;
				const auto current = queue.front();
				queue.pop();
				visited[current] = true;
				++nb_visited;

				for(unsigned b = 0; b < mer_ops::alpha; ++b) {
					const auto nmer = mer_ops::nmer(current, b);
					if(in_set(nmer)) continue; // Ignore mers in_set
					if(!visited[nmer]) {
						queue.push(nmer);
					} else if(nmer == m || nmer == rcm) {
						return true; // Loop involving rcm
					}
				}
			}

			if(visited[rcm]) break;
			queue.push(rcm);
		}

		return false;
	}

	void clear_queue() {
		while(!queue.empty()) queue.pop();
	}
};

// Greedy optimization procedure for mer_set using the random order. If
// can_super is true, result is a super set of the canonicalized mer_set.
struct greedy_opt {
	// tarjan_scc<mer_ops> comp_scc;
	symm_bfs bfs;
	std::unordered_set<mer_t> remove;

	template<typename C>
	std::pair<mer_t, mer_t> optimize(const C& mer_set, const std::vector<mer_t>& order, bool can_super, uint64_t max_iteration = 0) {
		remove.clear();
		// std::cout << "counts_orig" << std::endl;
		// const auto counts_orig = comp_scc.scc_counts(is_in_set(mer_set));

		uint64_t iteration = 0;
		for(const auto m : order) {
			if(can_super && m <= mer_ops::reverse_comp(m)) continue; // Must it be a super-set of canonical

			// Try adding m to remove. If increase SCCs, don't keep it
			remove.insert(m);
			const is_in_opt opt(mer_set, remove);
			std::cout << "remove " << (size_t)m << std::endl;
			// const auto counts = comp_scc.scc_counts(opt);
			// if(counts.first > counts_orig.first || counts.second > counts_orig.second)
			// 	remove.erase(m);
			if(bfs.has_cycle(opt, m))
				remove.erase(m);

			++iteration;
			std::cout << "iteration " << iteration << std::endl;
			if(max_iteration > 0 && iteration >= max_iteration) break;
		}

		return std::make_pair(0, 0);
	}
};

int main(int argc, char* argv[]) {
	opt_canon args(argc, argv);

    auto prg = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
                                           args.iseed_given ? args.iseed_arg : nullptr);

	std::cout << "Read set" << std::endl;
	const auto mer_set = get_mds<std::unordered_set<mer_t>>(args.sketch_file_arg, args.sketch_arg);
	std::cout << "Shuffle set" << std::endl;
	std::vector<mer_t> order(mer_set.begin(), mer_set.end());
	std::shuffle(order.begin(), order.end(), prg);

	greedy_opt optimizer;

	const auto counts = optimizer.optimize(mer_set, order, true, args.iteration_given? args.iteration_arg : 0);
	std::cout << "set\tsize\tsccs\n"
			  << "orig\t" << mer_set.size() << '\t' << counts << '\n'
			  << "union\t" << union_size(mer_set, (std::set<mer_t>*)nullptr) /* << '\t' << optimizer.comp_scc.scc_counts(is_in_union(mer_set)) */ << '\n'
			  << "canon\t" << canonicalize_size(mer_set) /* << '\t' << optimizer.comp_scc.scc_counts(can_is_in_set(mer_set)) */ << '\n'
			  << "super\t" << union_size(mer_set, &optimizer.remove) /* << '\t' << optimizer.comp_scc.scc_counts(is_in_opt(mer_set, optimizer.remove)) */ << '\n';

	optimizer.optimize(mer_set, order, false, args.iteration_given? args.iteration_arg : 0);
	std::cout << "opt\t" << union_size(mer_set, &optimizer.remove) /* << '\t' <<  optimizer.comp_scc.scc_counts(is_in_opt(mer_set, optimizer.remove)) */ << '\n';

	return EXIT_SUCCESS;
}
