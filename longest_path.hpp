#ifndef LONGEST_PATH_H_
#define LONGEST_PATH_H_

#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

#include "common.hpp"
#include "dbg.hpp"
#include "mds_op.hpp"

// Returns the longest path in number of k-mers/nodes in the de Bruijn graph
// decycled by bmds. (NOTE: bmds is copied as a vector.) So the longest string
// is longest + k - 1
template<typename mer_op_type>
struct longest_path_type {
    typedef typename mer_op_type::mer_t mer_t;
    std::vector<tristate_t>& visited;
    longest_path_type()
        : visited(mer_op_type::nb_mers, no)
        {}

    mer_t longest_path(const std::vector<tristate_t>& bmds, const std::vector<mer_t>& fms) {
        assert2(bmds.size() == visited.size(), "Wrong bmds input bmds size");
        std::copy(bmds.begin(), bmds.end(), visited.begin());
        std::stack<mer_t> indeg0;
        std::vector<mer_t> paths(mer_op_type::nb_mers, 1); // 1 for itself as a k-mer.
        mer_t longest = 0;

        for(const auto fm : fms) {
            for(mer_t b = 0; b < mer_op_type::nb_mers; ++b) {
                assert2(visited[mer_op_type::lc(fm, b)] == yes, "Invalid F-move, missing left companion " << fm << ' ' << b);
                assert2(visited[mer_op_type::rc(fm, b)] != yes, "Invalid F-move, right companion is present " << fm << ' ' << b);
                indeg0.push(mer_op_type::rc(fm, b));
            }
        }

        while(!indeg0.empty()) {
            const auto m = indeg0.top();
            indeg0.pop();
            visited[m] = yes;
            const bool push = mer_op_type::has_fm(visited, mer_op_type::fmove(m));

            for(mer_t b = 0; b < mer_op_type::nb_mers; ++b) {
                const auto nm = mer_op_type::nmer(m, b);
                if(visited[nm] == yes) continue;
                const auto len = std::max(paths[nm], paths[m] + 1);
                paths[nm] = len;
                longest = std::max(longest, len);
                if(push) indeg0.push(nm);
            }
        }

        return longest;
    }

    mer_t longest_path(const std::vector<mer_t>& mds) {
        std::vector<tristate_t> bmds;
        std::vector<mer_t> fms;

        mds_op_type<mer_op_type>::from_mds_fms(mds, bmds, fms);
        return longest_path_type(bmds, fms);
    }
};
#endif // LONGEST_PATH_H_
