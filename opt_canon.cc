#include <iostream>
#include <unistd.h>
#include <vector>
#include <unordered_set>
#include <random>
#include <chrono>
#include <bitset>
#include <cstdlib>
#include <csignal>
#include <atomic>
#include <functional>

#include "opt_canon.hpp"
#include "misc.hpp"
#include "common.hpp"
#include "mt_queue.hpp"
#include "random_seed.hpp"
#include "simple_thread_pool.hpp"


#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
typedef amer_type<K, ALPHA> amer_t;
typedef amer_t::mer_ops mer_ops;
typedef amer_t::mer_t mer_t;

// Does a BFS to detect a new cycle in the de Bruijn graph minus a set. Starts
// from m and the reverse complement of m (rcm) and check for a loop back to m
// or back to rcm.
template<typename mer_ops>
struct symm_bfs {
	mt_queue<amer_t> _queue;
	// std::vector<std::atomic<char>> _visited;
	std::vector<char> _visited; // Don't use atomic operations. See mark_visited().
	simple_thread_pool<std::function<void(int)>> _pool;

	static void noprogress(amer_t i) {}

	symm_bfs(int ths)
		: _queue(mer_ops::nb_mers)
		, _visited(mer_ops::nb_mers)
		, _pool(ths)
		{}
	~symm_bfs() { _pool.stop(); }

	template<typename Fn>
	bool has_cycle(Fn in_set, amer_t m) {
		// std::fill(_visited.begin(), _visited.end(), 0);
		std::memset(_visited.data(), 0, _visited.size() * sizeof(typename decltype(_visited)::value_type));
		_queue.clear();
		volatile bool found_loop = false;
		// The reverse complement of m is also considered removed from set and a
		// loop involving rcm also triggers returning true.
		const auto rcm = m.reverse_comp();

		// Process one when starting from m. Consider rcm not part of the set.n
		auto process_level = [&](int) {
			std::pair<amer_t*,ssize_t> push_loc{nullptr, 0};
			ssize_t push_index = 0;

			while(!found_loop) {
				const auto slice = _queue.multi_pop();
				if(slice.second <= 0) break; // Finished queue of current level

				// Slice of length slice.second or ends with sentinel value m
				for(ssize_t i = 0; i < slice.second && slice.first[i] != m; ++i) {
					amer_t::mer_rc_pair nmer_rc(slice.first[i].nmer(0));
					for(unsigned b = 0; b < mer_ops::alpha; ++b, ++nmer_rc) {
						if(mark_visited(nmer_rc.mer)) {
							const bool is_in_set = (nmer_rc.mer != rcm) && in_set(nmer_rc);
							if(is_in_set) continue; // Ignore mers in_set
							if(push_index >= push_loc.second) {
								push_loc = _queue.multi_push();
								push_index = 0;
							}
							push_loc.first[push_index] = nmer_rc.mer;
							++push_index;
						} else if(nmer_rc.mer == m) {
							found_loop = true; // Loop involving m or rcm
							break;
						}
					}
				}
			}

			// Padd unfilled location with sentinel value m
			if(push_loc.first && push_index < push_loc.second)
			    push_loc.first[push_index] = m;
		};
		_pool.set_work(process_level);

		// Prime queue. Simpler but equivalent to process_level
		mark_visited(m);
		amer_t::mer_rc_pair nmer_rc(m.nmer(0));
		for(unsigned b = 0; b < mer_ops::alpha; ++b, ++nmer_rc) {
			if(mark_visited(nmer_rc.mer)) {
				if((nmer_rc.mer != rcm) && in_set(nmer_rc)) continue;
				_queue.push(nmer_rc.mer);
			} else if(nmer_rc.mer == m) {
				return true;
			}
		}
		_queue.swap();

		while(!_queue.current_empty() && !found_loop) {
			_pool.start();
			_queue.swap();
		}

		return found_loop;
	}

	// Mark node m as _visited. Returns true if not previously visited. I.e.,
	// this call is the one who changed it to visited.
	bool mark_visited(const amer_t& m) {
		// Don't use any atomic operations to save time, although this is not
		// strictly correct for a BFS. Meaning a node could be visited multiple
		// times. This is rare. More importantly, it may add a bit of useless
		// work but it doesn't affect the correctness. Overall it is worth it.
        const auto prev = _visited[m.val];
		_visited[m.val] = 1;
		return prev == 0;
		// return _visited[m.val].exchange(1) == 0;
	}
};


struct is_in_set {
	const std::unordered_set<amer_t>& set;
	is_in_set(const std::unordered_set<amer_t>& s) : set(s) {}
	bool operator()(amer_t m) const { return set.find(m) != set.cend(); }
};

template<typename S>
struct is_in_union {
	const S& set;
	is_in_union(const S& s) : set(s) {}
	bool operator()(amer_t m) const {
		return set.find(m) != set.cend() || set.find(m.reverse_comp()) != set.cend();
	}
	bool operator()(const amer_t::mer_rc_pair& pair) const {
		return set.find(pair.mer) || set.find(pair.rc);
	}
};

template<typename C>
size_t canonicalize_size(const C& mers) {
	size_t size = 0;
	for(const auto& m : mers) {
		const auto rcm = m.reverse_comp();
		// Add two for canonical k-mers (1 for itself, 1 for its rc), unless it is self rc (then add only 1)
		if(m < rcm || m == rcm) {
			++size;
			if(m < rcm)
				++size;
		}
	}

	return size;
}

template<typename C, typename S>
size_t union_size(const C& mers, const S& set) {
	size_t size = 0;
	for(const auto& m : mers) {
		if(set.find(m) == set.end()) continue;
		++size; // Add one for the k-mer itself
		const auto rcm = m.reverse_comp();
		if(m != rcm && set.find(rcm) == set.end())
			++size; // Add one for its rc if not in set
	}
	return size;
}

// set as bitset for quick membership
template<typename mer_ops>
struct quickset {
	typedef amer_t value_type;
	std::bitset<mer_ops::nb_mers>* _data;
	quickset()
	: _data(new std::bitset<mer_ops::nb_mers>)
		{}
	~quickset() { delete _data; }

	void set(const amer_t& x) { _data->set(x.val); }
	void erase(const amer_t& x) { _data->reset(x.val); }

	bool find(const amer_t& x) const { return _data->test(x.val); }
	constexpr bool end() const { return false; }
	constexpr bool cend() const { return false; }
};

namespace
{
    volatile std::sig_atomic_t terminate = 0;
}
 
void signal_handler(int signal)
{
    terminate = 1;
}

template<typename mer_ops, bool enabled>
struct amain {
    int operator()(const opt_canon& args) {
        std::cerr << "Problem size too big" << std::endl;
		return EXIT_FAILURE;
    }
};

template<typename mer_ops>
struct amain<mer_ops, true> {
    int operator()(const opt_canon& args) {
		auto prg = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
											   args.iseed_given ? args.iseed_arg : nullptr);

		// Install a signal handler so the computation can be stopped at any time
		std::signal(SIGINT, signal_handler);
		std::signal(SIGTERM, signal_handler);



		auto orig_set = get_mds<std::unordered_set<amer_t>>(args.sketch_file_arg, args.sketch_arg);
		//auto order = get_mds<std::vector<amer_t>>(args.sketch_file_arg, args.sketch_arg);
		std::vector<amer_t> order(orig_set.cbegin(), orig_set.cend());
		std::shuffle(order.begin(), order.end(), prg);

		quickset<mer_ops> mer_set;
		for(const auto& m : order)
			mer_set.set(m);

		is_in_union union_set(mer_set);
		size_t removed = 0;

		int nb_threads = args.threads_arg > std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : args.threads_arg;
 		symm_bfs<mer_ops> bfs(nb_threads);
		std::cout << "original set: " << order.size()
				  << "\ncanonicalized set: " << canonicalize_size(order)
				  << "\nunion set: " << union_size(order, mer_set) << '\n';
		const auto begin = std::chrono::steady_clock::now();

		size_t progress = 0;
		const auto progress_suffix = isatty(1) ? '\r' : '\n';

		for(const auto& m : order) {
			if(terminate) break;
			if(args.progress_flag) {
				std::cout << progress << ' ' << removed << ' '
						  << (progress / (1e-6 + std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - begin).count()))
						  << progress_suffix << progress_suffix << std::flush;
				++progress;
			}

			const auto rcm = m.reverse_comp();
			// Must be a super set of canonicalize
			if(args.can_flag && (m < rcm || m == rcm)) continue;
			// Skip if rcm is also in set and not the canonical k-mer (avoid double computation)
			if(mer_set.find(rcm) != mer_set.end() && rcm < m) continue;
			const bool has_cycle = bfs.has_cycle(union_set, m) || bfs.has_cycle(union_set, rcm);
			if(!has_cycle) {
				++removed;
				mer_set.erase(m);
				mer_set.erase(rcm);
			}
		}
		if(progress) std::cout << '\n';
		std::cout << "Removed " << removed << '/' << progress << "\nunion set: " << union_size(order, mer_set) << '\n';

		if(args.output_given) {
			std::ofstream out(args.output_arg);
			bool first = true;
			for(mer_t i = 0; i < mer_ops::nb_mers; ++i) {
				if(!mer_set._data->test(i)) continue;
				if(!first) {
					out << ',';
					first = false;
				}
				out << amer_t(i);
				if(!out.good()) break;
			}
			out.close();
			if(!out.good()) {
				std::cerr << "Error while writing set to '" << args.output_arg << "''" << std::endl;
				return EXIT_FAILURE;
			}
		}

		return EXIT_SUCCESS;
    }
};



int main(int argc, char* argv[]) {
	opt_canon args(argc, argv);

return amain<mer_ops, mer_ops::ak_bits <= mer_ops::max_bits>()(args);
}
