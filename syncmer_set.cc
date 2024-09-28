#include <random>
#include "argparse.hpp"
#include "random_seed.hpp"
#include "divisor.hpp"
#include "syncmer.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

struct SyncmerArgs : argparse::Args {
	uint32_t& s_arg = kwarg("s", "s-mer size");
	uint32_t& t_arg = kwarg("t", "Position of minimum s-mer");
	std::optional<const char*>& iseed_arg = kwarg("i,iseed", "Input seed file");
	std::optional<const char*>& oseed_arg = kwarg("o,ioeed", "Output seed file");

	void welcome() override {
		std::cout << "Generate a syncmer set" << std::endl;
	}
};

int main(int argc, char* argv[]) {
	const auto args = argparse::parse<SyncmerArgs>(argc, argv);

	const auto nb_smers = ipow(mer_ops::alpha, args.s_arg);
	jflib::divisor64 div_s(nb_smers); // Fast division / modulo by nb_smers

	// The properly seeded PRG
	auto prg = seeded_prg<std::mt19937_64>(args.oseed_arg ? *args.oseed_arg : nullptr,
                                           args.iseed_arg ? *args.iseed_arg : nullptr);

	// Random order of s-mers
	std::vector<mer_t> smer_order(nb_smers);
	for(mer_t s = 0; s < nb_smers; ++s)
		smer_order[s] = s;
	std::shuffle(smer_order.begin(), smer_order.end(), prg);

	bool first = true;
	for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
		if(min_smer<mer_ops>(m, args.s_arg, smer_order, div_s) == args.t_arg) {
			if(!first) {
				std::cout << ',';
			} else {
				first = false;
			}
			std::cout << (size_t)m;
		}
	}
	std::cout << '\n';

	return EXIT_SUCCESS;
}
