#include <random>
#include "syncmer_set.hpp"
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



int main(int argc, char* argv[]) {
	syncmer_set args(argc, argv);
	const auto nb_smers = ipow(mer_ops::alpha, args.s_arg);
	jflib::divisor64 div_s(nb_smers); // Fast division / modulo by nb_smers

	// The properly seeded PRG
	auto prg = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
                                           args.iseed_given ? args.iseed_arg : nullptr);

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
