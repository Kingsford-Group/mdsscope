#include <cstdlib>
#include <unordered_set>

#include "misc.hpp"
#include "sequence.hpp"
#include "sketch_histo.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

typedef bool (*lookup_fn)(const std::unordered_set<mer_t>& set, mer_t m);

bool straight_set(const std::unordered_set<mer_t>& set, mer_t m) {
	return set.find(m) != set.end();
}

bool canonical_set(const std::unordered_set<mer_t>& set, mer_t m) {
	return set.find(mer_ops::canonical(m)) != set.end();
}

bool union_set(const std::unordered_set<mer_t>& set, mer_t m) {
	return (set.find(m) != set.end()) || (set.find(m) != set.end());
}

int main(int argc, char* argv[]) {
	sketch_histo args(argc, argv);
	const auto mer_set = get_mds<std::unordered_set<mer_t>>(args.sketch_file_arg, args.sketch_arg);
	std::cout << "mer_set size " << mer_set.size() << '\n';

	char inchar = '0';
	size_t offset = 0;
	mer_t mer = 0;
	translated_stream ts(args.alphabet_arg, mer_ops::alpha, std::cin);
	// Read first s-1 bases
	for( ; offset + 1 < mer_ops::k && ts >> inchar; ++offset)
		mer = mer_ops::nmer(mer, inchar);

	lookup_fn lookup = straight_set;
	if(args.canonical_flag)
		lookup = canonical_set;
	else if(args.union_flag)
		lookup = union_set;

	std::vector<size_t> histo;
	size_t prev = 0;
	offset = 0;
	while(ts >> inchar) {
		mer = mer_ops::nmer(mer, inchar);
		if(lookup(mer_set, mer)) {
			// std::cout << ' ' << offset << ' ' << prev;
			const size_t dist = offset - prev;
			if(dist >= histo.size())
				histo.resize(dist + 1, 0);
			++histo[dist];
			prev = offset;
		}
		// std::cout << '\n';
		++offset;
	}
	if(!ts.is.eof()) {
		std::cerr << "Encountered error reading sequence" << std::endl;
		return EXIT_FAILURE;
	}

	for(size_t i = 0; i < histo.size(); ++i) {
		if(histo[i])
			std::cout << i << ' ' << histo[i] << '\n';
	}

	return EXIT_SUCCESS;
}