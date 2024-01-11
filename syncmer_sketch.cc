#include <iostream>
#include <random>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "random_seed.hpp"
#include "syncmer_sketch.hpp"
#include "divisor.hpp"
#include "sequence.hpp"
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
	syncmer_sketch args(argc, argv);
	if(args.s_arg > mer_ops::k) {
		std::cerr << "Error: s must be less than k" << std::endl;
		return EXIT_FAILURE;
	}
	// const unsigned nb_s_mers = mer_ops::k - args.s_arg + 1;
	const auto nb_smers = ipow(mer_ops::alpha, args.s_arg);
	jflib::divisor64 div_s(nb_smers); // Fast division / modulo by nb_smers

	// The properly seeded PRG
	auto prg = seeded_prg<std::mt19937_64>(args.oseed_given ? args.oseed_arg : nullptr,
                                           args.iseed_given ? args.iseed_arg : nullptr);
	// Random order of s-mers
	std::vector<mer_t> smer_order(nb_smers);
	for(mer_t s = 0; s < nb_smers; ++s)
		smer_order[s] = s;
	if(!args.lex_flag)
		std::shuffle(smer_order.begin(), smer_order.end(), prg);

	char inchar = '0';
	mer_t mer = 0;
	size_t offset = 0;
	translated_stream ts(args.alphabet_arg, mer_ops::k, std::cin);
	// Read first s-1 bases
	while(offset + 1 < mer_ops::k && ts >> inchar) {
		mer = mer_ops::nmer(mer, inchar);
		// std::cout << "s-mer " << (size_t)inchar << ' ' << (size_t)mer << '\n';
		++offset;
	}

	std::vector<size_t> histo;
	size_t prev = 0;
	offset = 0;
	while(ts >> inchar) {
		mer = mer_ops::nmer(mer, inchar);
		const unsigned min = min_smer<mer_ops>(args.canonical_flag ? mer_ops::canonical(mer) : mer, args.s_arg, smer_order, div_s);
		if(min == args.t_arg) {
			// std::cout << (size_t)mer << ' ' << offset << ' ' << prev << '\n';
			const size_t dist = offset - prev;
			if(dist >= histo.size())
				histo.resize(dist + 1, 0);
			++histo[dist];
			prev = offset;
		}
		++offset;
	}

	for(size_t i = 0; i < histo.size(); ++i) {
		if(histo[i])
			std::cout << i << ' ' << histo[i] << '\n';
	}

	return EXIT_SUCCESS;
}
