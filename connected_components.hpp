#ifndef CONNECTED_COMPONENTS_H_
#define CONNECTED_COMPONENTS_H_

#include "common.hpp"
#include "mds_op.hpp"
#include "union_rank.hpp"

// Returns the number of connected components after decycling by bmds. The
// assignment to components is in the component union rank data structure.
//
// XXX: too much code overlap with longest_path. Refactor
template<typename mer_op_type>
struct connected_components_type {
	typedef typename mer_op_type::mer_t mer_t;
	typedef mds_op_type<mer_op_type> mds_ops;

	union_rank<mer_op_type> components;

	connected_components_type() {}

	void compute_components(const std::vector<tristate_t>& bmds) {
		for(mer_t m = 0; m < mer_op_type::nb_mers; ++m) {
			// std::cout << (size_t)m << ' ' << (int)bmds[m] << '\n';
			if(bmds[m] == yes) continue;

            for(mer_t b = 0; b < mer_op_type::alpha; ++b) {
                const auto nm = mer_op_type::nmer(m, b);
				if(bmds[nm] == yes) continue;
				// std::cerr << "merge " << (size_t)m << ' ' << (size_t)nm << '\n';
				components.merge(m, nm);
            }
        }
	}

	mer_t nb_components(const std::vector<mer_t>& mds) {
		std::vector<tristate_t> bmds;
		std::vector<mer_t> fms;

		mds_op_type<mer_op_type>::from_mds_fms(mds, bmds, fms);
		compute_components(bmds);
		return components.nb_sets - mds.size();
	}
};

#endif // CONNECTED_COMPONENTS_H_
