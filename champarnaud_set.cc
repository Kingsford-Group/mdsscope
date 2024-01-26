#include <cstdlib>
#include <vector>
#include <bitset>

#include "divisor.hpp"
#include "misc.hpp"
#include "common.hpp"
#include "dbg.hpp"
#include "champarnaud.hpp"
#include "champarnaud_set.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

template<typename Fn>
void enumerate_champarnaud(Fn fn) {
    index_t<mer_ops> index(K);

	std::vector<bool> visited(mer_ops::nb_mers);
	for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
		if(visited[m]) continue;
		for(mer_t cur = m; !visited[cur]; cur = mer_ops::nmer(cur)) // Mark PCR as visited
			visited[cur] = true;
		fn(index.reversed(m));
	}

}

template<typename Fn>
void brute_champarnaud(Fn fn) {
	index_t<mer_ops> index(K);

	for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
		if(index.in_champarnaud_set(m))
			fn(m);
	}
}

int main(int argc, char* argv[]) {
	champarnaud_set args(argc, argv);

    if(args.notsorted_flag) {
        bool first = true;
        auto fn = [&first](mer_t m) {
            if(!first) {
                std::cout << ',';
            } else {
                first = false;
            }
            std::cout << (size_t)m;
        };
		!args.brute_flag ? enumerate_champarnaud(fn) : brute_champarnaud(fn);
    } else {
        std::vector<mer_t> set;
        auto fn = [&set](mer_t m) { set.push_back(m); };
		!args.brute_flag ? enumerate_champarnaud(fn) : brute_champarnaud(fn);
        std::sort(set.begin(), set.end());
        std::cout << joinT<size_t>(set, ',');
    }

    std::cout << '\n';

	return EXIT_SUCCESS;
}
