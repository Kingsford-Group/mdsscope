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

struct alphabet_t {
	std::vector<char> table;
	alphabet_t(bool custom, const char* str) : table(256, -1) {
		if(custom) {
			if(strlen(str) != mer_ops::alpha)
				throw std::runtime_error("Alphabet translation of wrong length");
			for(unsigned i = 0; i < strlen(str); ++i) {
				if(!std::isprint(str[i]) || std::isspace(str[i]))
					throw std::runtime_error("Invalid character in alphabet");
				table[str[i]] = i;
			}
		} else {
			for(unsigned i = 0; i < mer_ops::alpha; ++i)
				table['0' + i] = i;
		}
	}

	char operator[](char c) const {
		char res = table[c];
		if(res < 0) [[unlikely]]
			throw std::runtime_error(std::string("Invalid character in input: ") + res);
	return res;
	}


};

struct translated_stream {
	const alphabet_t& alphabet;
	std::istream& is;
	size_t offset;
	translated_stream(const alphabet_t& a, std::istream& i)
		: alphabet(a)
		, is(i)
		, offset(0)
		{}

	operator bool() const { return (bool)is; }

	translated_stream& operator>>(char& c) {
		char rc;
		while(is >> rc) {
			++offset;
			// std::cout << "read " << offset << ' ' << rc << ' ' << std::isspace(rc) << ' ' << (int)alphabet.table[rc] << '\n';
			if(std::isspace(rc)) continue;
			c = alphabet.table[rc];
			if(c < 0) [[unlikely]] {
				std::ostringstream msg;
				msg << "Invalid character '" << rc << "' at position " << offset;
				throw std::runtime_error(msg.str());
			}
			break;
		}
		return *this;
	}
};

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

	alphabet_t alphabet(args.alphabet_given, args.alphabet_arg);

	char inchar = '0';
	mer_t mer = 0;
	size_t offset = 0;
	translated_stream ts(alphabet, std::cin);
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
