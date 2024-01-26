#ifndef SYNCMER_H_
#define SYNCMER_H_

#include<vector>
#include "divisor.hpp"

template<typename mer_ops, typename mer_t = typename mer_ops::mer_t>
unsigned min_smer(mer_t m, unsigned s, const std::vector<mer_t>& smer_order, const jflib::divisor64& div_s) {
	// std::cout << "min_smer " << (size_t)m;
	unsigned min = 0;
	mer_t min_val = smer_order[m % div_s];
	m /= mer_ops::alpha;
	for(unsigned i = 1; i < mer_ops::k - s + 1; ++i, m /= mer_ops::alpha) {
		mer_t smer = smer_order[m % div_s];
		if(smer <= min_val) {
			min = i;
			min_val = smer;
		}
	}
	// std::cout << ' ' << mer_ops::k - s - min << '\n';
	return mer_ops::k - s - min;
}


#endif // SYNCMER_H_
