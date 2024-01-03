#include <vector>
#include <random>
#include <numeric>

#include "random_seed.hpp"
#include "common.hpp"
#include "frac_set.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

constexpr mer_t nb_necklace() {
	mer_t res = 0;
	for(unsigned i = 1; i <= mer_ops::k; ++i) {
		res += ipow(mer_ops::alpha, std::gcd(i, mer_ops::k));
	}
	return res / mer_ops::k;
}

int main(int argc, char* argv[]) {
	frac_set args(argc, argv);

	size_t set_size;
	if(args.fraction_given) {
		if(args.fraction_arg < 0.0 || args.fraction_arg > 1.0) {
			std::cerr << "Frace must be in [0, 1]" << std::endl;
			return EXIT_FAILURE;
		}
		set_size = std::min((size_t)mer_ops::nb_mers, (size_t)std::round(args.fraction_arg * mer_ops::nb_mers));
	} else if(args.size_given) {
		if(args.size_arg > mer_ops::nb_mers) {
			std::cerr << "Size must be in [0, " << (size_t)mer_ops::nb_mers << ']' << std::endl;
			return EXIT_FAILURE;
		}
		set_size = args.size_arg;
	} else {
		set_size = nb_necklace();
	}
	std::cerr << "size " << set_size << '\n';

	// The properly seeded PRG
	auto prg = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
                                           args.iseed_given ? args.iseed_arg : nullptr);

	std::vector<mer_t> all_mers(mer_ops::nb_mers);
	for(mer_t m = 0; m < all_mers.size(); ++m)
		all_mers[m] = m;
	std::shuffle(all_mers.begin(), all_mers.end(), prg);

	std::cout << joinitT<size_t>(all_mers.begin(), all_mers.begin() + set_size, ',') << '\n';

	return EXIT_SUCCESS;
}
