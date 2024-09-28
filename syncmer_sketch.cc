#include <iostream>
#include <random>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include "argparse.hpp"
#include "random_seed.hpp"
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

struct SyncmerSketchArgs : argparse::Args {
	uint32_t& s_arg = kwarg("s", "s-mer size");
	uint32_t& t_arg = kwarg("t", "Position of minimum s-mer");
	std::optional<const char*>& alphabet_arg = kwarg("a,alphabet", "Alphabet translation");
	bool& canonical_flag = flag("c,canonical", "Canonical space");
	bool& lex_flag = flag("lex", "Use lexicographic order");
	std::optional<const char*>& iseed_arg = kwarg("i,iseed", "Input seed file");
	std::optional<const char*>& oseed_arg = kwarg("o,ioeed", "Output seed file");
};

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;


int main(int argc, char* argv[]) {
	const auto args = argparse::parse<SyncmerSketchArgs>(argc, argv);

	if(args.s_arg > mer_ops::k) {
		std::cerr << "Error: s must be less than k" << std::endl;
		return EXIT_FAILURE;
	}
	// const unsigned nb_s_mers = mer_ops::k - args.s_arg + 1;
	const auto nb_smers = ipow(mer_ops::alpha, args.s_arg);
	jflib::divisor64 div_s(nb_smers); // Fast division / modulo by nb_smers

	// The properly seeded PRG
	auto prg = seeded_prg<std::mt19937_64>(args.oseed_arg ? *args.oseed_arg : nullptr,
                                           args.iseed_arg ? *args.iseed_arg : nullptr);
	// Random order of s-mers
	std::vector<mer_t> smer_order(nb_smers);
	for(mer_t s = 0; s < nb_smers; ++s)
		smer_order[s] = s;
	if(!args.lex_flag)
		std::shuffle(smer_order.begin(), smer_order.end(), prg);

	char inchar = '0';
	mer_t mer = 0;
	size_t offset = 0;
	translated_stream ts(*args.alphabet_arg, mer_ops::k, std::cin);
	// Read first k-1 bases
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
