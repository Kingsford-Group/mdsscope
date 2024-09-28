#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <vector>
#include <algorithm>
#include <random>
#include <memory>

#include "argparse.hpp"
#include "misc.hpp"
// #include "common.hpp"
#include "permutations.hpp"
#include "random_seed.hpp"
#include "sequence.hpp"
#include "mykkeltveit.hpp"
#include "champarnaud.hpp"
#include "syncmer.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

struct SketchHistoArgs : argparse::Args {
	std::optional<const char*>& alphabet_arg = kwarg("a,alphabet", "Alphabet translation");
	std::optional<const char*>& sketch_file_arg = kwarg("f,sketch-file", "File with sketch mer set");
	std::optional<const char*>& iseed_arg = kwarg("i,iseed", "Input seed file");
	std::optional<const char*>& oseed_arg = kwarg("o,ioeed", "Output seed file");
	bool& mykkeltveit_flag = flag("mykkeltveit", "Stream Mykkeltveit set");
	std::optional<uint32_t>& syncmer_arg = kwarg("syncmer", "Stream syncmer set");
	std::optional<uint32_t>& syncmer_s_arg = kwarg("syncmer-s", "Syncmer s parameter (default k/2 - 1)");
	std::optional<double>& frac_arg = kwarg("frac", "Stream Frac set");
	bool& champarnaud_flag = flag("champarnaud", "Stream Champarnaud set");
	bool& straight_flag = flag("s,straight", "Use set directly (default)");
	bool& canonical_flag = flag("c,canonical", "Use canonical k-mers");
	bool& union_flag = flag("u,union", "Use union of set and reverse complement set");
	bool& sum_flag = flag("sum", "Output weighted sum of the histo");
	uint32_t& hmin_arg = flag("hmin", "Ignore length of histo less than hmin").set_default(0);


	std::vector<const char*>& sketch_arg = arg("sketch").set_default("");

	void welcome() override {
		std::cout <<
			"Sketch a sequence given a context free / set scheme\n\n"
			"Boils down to the intersection of k-mers in the set and in the sequence.\n"
			"Sequence is read from stdin, set from -f or command line arguments."
			<< std::endl;
	}
};

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;


bool straight_set(const std::unordered_set<mer_t>* set, mer_t m) {
	return set->find(m) != set->end();
}

// struct memoized {
// 	std::unordered_map<mer_t,bool> cache;
// 	std::function<bool(mer_t)> f;
// 	memoized(std::function<bool(mer_t)> f) : f(f) {}
// 	bool operator()(mer_t m) {
// 		auto it = cache.find(m);
// 		if(it != cache.end())
// 			return it->second;
// 		const bool v = f(m);
// 		cache.insert(it, std::make_pair(m, v));
// 		return v;
// 	}
// };

bool memoized(std::unordered_map<mer_t,bool>* cache, std::function<bool(mer_t)> f, mer_t m) {
	auto it = cache->find(m);
	if(it != cache->end())
		return it->second;
	const bool v = f(m);
	cache->insert(it, std::make_pair(m, v));
	return v;
}

bool mykkeltveit(const root_unity_type<mer_ops>* root_unity, mer_t m) {
	return root_unity->in_mykkeltveit_set(m);
}
// struct mykkeltveit {
// 	const root_unity_type<mer_ops> root_unity;

// 	inline bool operator()(mer_t m) { return root_unity.in_mykkeltveit_set(m); }
// };

struct champarnaud_data_type {
	const index_t<mer_ops> index;
	champarnaud_data_type() : index(mer_ops::k) {}
};
inline bool champarnaud(const champarnaud_data_type* d, mer_t m) { return d->index.in_champarnaud_set(m); }

// Syncmer when s is small (order on s-mer is an explicitely shuffled array)
struct syncmer_data_type {
	const unsigned s, t;
	const jflib::divisor64 div_s;
	std::vector<mer_t> smer_order;

	template<typename PRG>
	syncmer_data_type(unsigned s, unsigned t, PRG* prg)
		: s(s)
		, t(t)
		, div_s(ipow((mer_t)mer_ops::alpha, (mer_t)s))
		, smer_order(div_s.d())
		{
			for(mer_t s = 0; s < div_s.d(); ++s)
				smer_order[s] = s;
			if(prg)
				std::shuffle(smer_order.begin(), smer_order.end(), *prg);

		}
};
inline bool syncmer(const syncmer_data_type* d, mer_t m) {
	return min_smer<mer_ops>(m, d->s, d->smer_order, d->div_s) == d->t;
}

// Syncmer for large s (s-mer order using perm)
struct syncmer_large_data_type {
    const unsigned s, t;
	const jflib::divisor64 div_s;
    const LubyRackofPermutation<mer_t> perm;

	template<typename PRG>
	syncmer_large_data_type(unsigned s, unsigned t, PRG* prg)
	: s(s)
	, t(t)
	, div_s(ipow((mer_t)mer_ops::alpha, (mer_t)s))
	, perm(*prg)
	{ }
};
inline bool syncmer_large(const syncmer_large_data_type* d, mer_t m) {
	return min_large_smer<mer_ops>(m, d->s, d->perm, d->div_s) == d->t;
}

// struct frac_data_type {
// 	typedef std::uniform_int_distribution<mer_t> mask_rng;
// 	const mer_t mask1, mask2;
// 	// jflib::divisor64 rot;
// 	mer_t rot;
// 	const mer_t shift;
// 	const mer_t thresh;
// 	template<typename PRG>
// 	frac_data_type(double f, PRG& prg)
// 		: mask1(mask_rng(0, std::numeric_limits<mer_t>::max())(prg))
// 		, mask2(mask_rng(0, std::numeric_limits<mer_t>::max())(prg))
// 		, rot(ipow((mer_t)mer_ops::alpha, std::uniform_int_distribution<int>(0, mer_ops::k-1)(prg)))
// 		// , shift(mer_ops::nb_mers / rot.d())
// 		, shift(mer_ops::nb_mers / rot)
// 		, thresh(std::round(mer_ops::nb_mers * f))
// 		{}
// };
// inline bool frac(const frac_data_type* d, mer_t m) {
// 	m ^= d->mask1;
// 	m = (m % d->rot) * d->shift + (m / d->rot);
// 	m ^= d->mask2;
// 	m %= mer_ops::nb_mers;
// 	return m < d->thresh;
// }

struct frac_data_type {
  const LubyRackofPermutation<mer_t> perm;
  const mer_t thresh;

  template<typename PRG>
  frac_data_type(double f, PRG& prg)
  : perm(prg)
  , thresh(std::round(std::pow(2.0, sizeof(mer_t) * 8) * f))
  { }
};

inline bool frac(const frac_data_type* d, mer_t m) {
  const auto val = d->perm(m);
  return  val < d->thresh;
}

bool canonical_fn(std::function<bool(mer_t)> f, mer_t m) {
	return f(mer_ops::canonical(m));
}

bool union_fn(std::function<bool(mer_t)> f, mer_t m) {
	return f(m) || f(mer_ops::canonical(m));
}

double fill_in_histo(translated_stream& ts, std::vector<size_t>& histo, std::function<bool(mer_t)> lookup) {
	std::fill(histo.begin(), histo.end(), 0);
	uint64_t selected = 0, kmers = 0;
	ts.header(); // Call at beginnin of every subsequence.
	// std::cout << "Reading " << ts.seq_name() << std::endl;

	mer_t mer = 0;
	char inchar = '0';
	size_t prev = 0;;
	{ // Read first s-1 bases
	  	size_t offset = 0;
		for( ; offset + 1 < mer_ops::k && ts >> inchar; ++offset) {
//			std::cout << "inchar " << (int)inchar << '\n';
			if(inchar == mer_ops::alpha) return -1.0;
			mer = mer_ops::nmer(mer, inchar);
		}
	}

	size_t offset = 0;
	while(ts >> inchar) {
		++kmers;
//		std::cout << "inchar " << (int)inchar << '\n';
		if(inchar == mer_ops::alpha) return -1.0;
		mer = mer_ops::nmer(mer, inchar);
//		std::cout << "lookup" << std::endl;
		if(lookup(mer)) {
			++selected;
//			std::cout << "true" << std::endl;
			const size_t dist = offset - prev;
			if(dist >= histo.size())
				histo.resize(dist + 1, 0);
			++histo[dist];
			prev = offset;
		}
//		std::cout << "after" << std::endl;
		++offset;
	}

	return (double)selected / (double)kmers;
}

int main(int argc, char* argv[]) {
	const auto args = argparse::parse<SketchHistoArgs>(argc, argv);

	std::vector<std::function<bool(mer_t)>> lookups; // Stack of lookup function. Composed with top of stack
	lookups.reserve(10);

	auto prg = seeded_prg<std::mt19937_64>(args.oseed_arg ? *args.oseed_arg : nullptr,
                                           args.iseed_arg ? *args.iseed_arg : nullptr);

	// Bottom layer
	std::unique_ptr<std::unordered_set<mer_t>> mer_set;
	std::unordered_map<mer_t,bool> mer_set_cache;
	std::unique_ptr<root_unity_type<mer_ops>> root_unity;
	std::unique_ptr<syncmer_data_type> syncmer_data;
	std::unique_ptr<syncmer_large_data_type> syncmer_large_data;
	std::unique_ptr<frac_data_type> frac_data;
	std::unique_ptr<champarnaud_data_type> champarnaud_data;

	if(args.sketch_file_arg || !args.sketch_arg.empty()) {
		mer_set.reset(new std::unordered_set<mer_t>);
		*mer_set = get_mds<std::unordered_set<mer_t>>(*args.sketch_file_arg, args.sketch_arg);
//	std::cout << "mer_set size " << mer_set.size() << ": " << joinT<size_t>(mer_set, ',') << '\n';
		lookups.emplace_back(std::bind_front(straight_set, mer_set.get())); // Bottom layer is querying the set
	} else if(args.mykkeltveit_flag) {
		root_unity.reset(new root_unity_type<mer_ops>);
		lookups.emplace_back(std::bind_front(mykkeltveit, root_unity.get()));
	} else if(args.syncmer_arg) {
		const unsigned s = args.syncmer_s_arg ? *args.syncmer_s_arg : K / 2 - 1;
		if(std::pow(mer_ops::alpha, s) < 1e9) {
            std::cerr << "syncmer " << s << ' ' << *args.syncmer_arg << std::endl;
            syncmer_data.reset(new syncmer_data_type(s, *args.syncmer_arg, &prg));
            lookups.emplace_back(std::bind_front(syncmer, syncmer_data.get()));
		} else {
            std::cerr << "syncmer large" << s << ' ' << *args.syncmer_arg << std::endl;
			syncmer_large_data.reset(new syncmer_large_data_type(s, *args.syncmer_arg, &prg));
			lookups.emplace_back(std::bind_front(syncmer_large, syncmer_large_data.get()));
		}
	} else if(args.frac_arg) {
		frac_data.reset(new frac_data_type(*args.frac_arg, prg));
		lookups.emplace_back(std::bind_front(frac, frac_data.get()));
	} else if(args.champarnaud_flag) {
		champarnaud_data.reset(new champarnaud_data_type());
		lookups.emplace_back(std::bind_front(champarnaud, champarnaud_data.get()));
	} else {
		std::cerr << "Missing set" << std::endl;
		return EXIT_FAILURE;
	}

	if(args.canonical_flag) {
	 	lookups.emplace_back(std::bind_front(canonical_fn, lookups.back()));
	} else if(args.union_flag) {
		lookups.emplace_back(std::bind_front(union_fn, lookups.back()));
	}

	lookups.emplace_back(std::bind_front(memoized, &mer_set_cache, lookups.back()));

	const auto& lookup = lookups.back();
	std::vector<size_t> histo;
	translated_stream ts(args.alphabet_arg ? *args.alphabet_arg : nullptr, mer_ops::alpha, std::cin);

	while(ts) {
		// std::cout << "loop" << std::endl;
		const auto density = fill_in_histo(ts, histo, lookup);
		if(!ts.seq_name().empty())
			std::cout << '>' << ts.seq_name() << '\n';
		std::cout << "# density " << density << '\n';
		if(args.sum_flag) {
			size_t sum = 0;
			for(size_t i = 0; i < histo.size(); ++i) {
				if(histo[i] && i >= args.hmin_arg)
					sum += i * histo[i];
			}
			std::cout << sum << '\n';
		} else {
			for(size_t i = 0; i < histo.size(); ++i) {
				if(histo[i] && i >= args.hmin_arg)
					std::cout << i << ' ' << histo[i] << '\n';
			}
		}
	}

	if(!ts.is.eof()) {
		std::cerr << "Encountered error reading sequence" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
