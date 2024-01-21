#ifndef TARJAN_SCC_H_
#define TARJAN_SCC_H_

#include <vector>
#include <utility>

template<typename mer_ops>
struct tarjan_scc {
	typedef typename mer_ops::mer_t mer_t;
	static constexpr mer_t undefined = mer_ops::nb_mers;
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


#endif // TARJAN_SCC_H_
