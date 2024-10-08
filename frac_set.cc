#include <vector>
#include <random>

#include "argparse.hpp"
#include "random_seed.hpp"
#include "common.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

struct FracSetArgs : argparse::Args {
	std::optional<double>& fraction_arg = kwarg("f,fraction", "Fraction of k-mers to keep");
	std::optional<double>& size_arg = kwarg("s,size", "Size of fractional set");
	std::optional<const char*>& iseed_arg = kwarg("i,iseed", "Input seed file");
	std::optional<const char*>& oseed_arg = kwarg("o,ioeed", "Output seed file");

};

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

int main(int argc, char* argv[]) {
	const auto args = argparse::parse<FracSetArgs>(argc, argv);

	size_t set_size;
	if(args.fraction_arg) {
		if(*args.fraction_arg < 0.0 || *args.fraction_arg > 1.0) {
			std::cerr << "Frace must be in [0, 1]" << std::endl;
			return EXIT_FAILURE;
		}
		set_size = std::min((size_t)mer_ops::nb_mers, (size_t)std::round(*args.fraction_arg * mer_ops::nb_mers));
	} else if(args.size_arg) {
		if(*args.size_arg > mer_ops::nb_mers) {
			std::cerr << "Size must be in [0, " << (size_t)mer_ops::nb_mers << ']' << std::endl;
			return EXIT_FAILURE;
		}
		set_size = *args.size_arg;
	} else {
		set_size = mer_ops::nb_necklaces;
	}
	std::cerr << "size " << set_size << '\n';

	// The properly seeded PRG
	auto prg = seeded_prg<std::mt19937_64>(args.oseed_arg ? *args.oseed_arg : nullptr,
                                           args.iseed_arg ? *args.iseed_arg : nullptr);

	std::vector<mer_t> all_mers(mer_ops::nb_mers);
	for(mer_t m = 0; m < all_mers.size(); ++m)
		all_mers[m] = m;
	std::shuffle(all_mers.begin(), all_mers.end(), prg);

	std::cout << joinitT<size_t>(all_mers.begin(), all_mers.begin() + set_size, ',') << '\n';

	return EXIT_SUCCESS;
}
