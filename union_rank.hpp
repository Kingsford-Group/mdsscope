#ifndef UNION_RANK_H_
#define UNION_RANK_H_

#include <vector>

#include "dbg.hpp"

// Disjoint sets / union-rank data structure

template<typename mer_op_type>
struct union_rank {
	typedef typename mer_op_type::mer_t mer_t;
	std::vector<mer_t> parent;
	std::vector<unsigned char> rank;
	mer_t nb_sets;

	union_rank()
		: parent(mer_op_type::nb_mers)
		, rank(mer_op_type::nb_mers, 0)
		, nb_sets(mer_op_type::nb_mers)
		{
			for(mer_t m = 0; m < mer_op_type::nb_mers; ++m)
				parent[m] = m; // Each mer alone in its set
		}

	mer_t find(mer_t x) {
		assert2(x < mer_op_type::nb_mers, "Invalid mer");

		while(parent[x] != x) {
			const auto tmp = parent[x];
			parent[x] = parent[parent[x]];
			x = tmp;
		}

		return x;
	}

	void merge(mer_t x, mer_t y) {
		x = find(x);
		y = find(y);

		if(x == y) return; // Already same set

		--nb_sets; // Actual merging
		if(rank[x] < rank[y])
			std::swap(x, y);

		parent[y] = x;
		if(rank[x] == rank[y])
			++rank[x];
	}
};

#endif // UNION_RANK_H_
