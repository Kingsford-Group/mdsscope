#ifndef LONGEST_PATH_H_
#define LONGEST_PATH_H_

#include <memory>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

#include "common.hpp"
#include "misc.hpp"
#include "mds_op.hpp"
#include "dbg.hpp"
#include "mds_op.hpp"

// Returns the longest path in number of k-mers/nodes in the de Bruijn graph
// decycled by bmds. (NOTE: bmds is copied as a vector.) So the longest string
// is longest + k - 1
template<typename mer_op_type>
struct longest_path_type {
    typedef typename mer_op_type::mer_t mer_t;
    typedef mds_op_type<mer_op_type> mds_ops;
    typedef std::set<std::vector<mer_t>> origin_set_type;
    typedef typename origin_set_type::iterator origin_it;

    std::vector<tristate_t> visited;
    std::vector<mer_t> paths;
    origin_set_type origin_set;
    std::vector<origin_it> origin;
    longest_path_type()
        : visited(mer_op_type::nb_mers, no)
        , paths(mer_op_type::nb_mers, 0)
        , origin(mer_op_type::nb_mers, origin_set.end())
        {}

    mer_t longest_path(const std::vector<tristate_t>& bmds, const std::vector<mer_t>& fms) {
        assert2(bmds.size() == visited.size(), "Wrong bmds input bmds size");
        std::copy(bmds.begin(), bmds.end(), visited.begin());
        origin_set.clear();

        std::vector<mer_t> indeg0;
        std::fill(paths.begin(), paths.end(), 1); // 1 for itself as a k-mer.
        mer_t longest = 1;

        for(const auto fm : fms) {
            for(mer_t b = 0; b < mer_op_type::alpha; ++b) {
                assert2(visited[mer_op_type::lc(fm, b)] == yes, "Invalid F-move, missing left companion " << fm << ' ' << b);
                const auto nm = mer_op_type::nmer(fm, b);
                if(visited[nm] != yes) {
                    indeg0.push_back(nm);
                    auto res = origin_set.insert(std::vector<mer_t>{nm});
                    origin[nm] = res.first;
                }
            }
        }

        while(!indeg0.empty()) {
            const auto m = indeg0.back();
            indeg0.pop_back();
            visited[m] = yes;
            const bool has_fm = mds_ops::has_fm(visited, mer_op_type::fmove(m));

            for(mer_t b = 0; b < mer_op_type::alpha; ++b) {
                const auto nm = mer_op_type::nmer(m, b);
                if(visited[nm] == yes) continue;
                const mer_t nlen = paths[m] + 1;
                if(nlen > paths[nm]) { // Strictly longer path. Origin via m
                    paths[nm] = nlen;
                    origin[nm] = origin[m];
                    longest = std::max(longest, nlen);
                } else if(nlen == paths[nm]) { // Same length path: new origin
                    std::vector<mer_t> norigin(*origin[nm]);
                    norigin.push_back(m);
                    std::sort(norigin.begin(), norigin.end());
                    auto res = origin_set.insert(norigin);
                    origin[nm] = res.first;
                }
                if(has_fm)
                    indeg0.push_back(nm);
            }
        }

        if(std::any_of(visited.begin(), visited.end(), [](tristate_t x) { return x != yes; }))
            return mer_op_type::nb_mers;

        return longest;
    }

    mer_t longest_path(const std::vector<mer_t>& mds) {
        std::vector<tristate_t> bmds;
        std::vector<mer_t> fms;

        mds_op_type<mer_op_type>::from_mds_fms(mds, bmds, fms);
        return longest_path(bmds, fms);
    }

    mer_t shortest_path(const std::vector<tristate_t>& bmds, const std::vector<mer_t>& fms) {
        // Because it is a DAG, do a BFS. Start from every F-move and traverse
        // the DAG until encountering right-companions (an RF-move). Mark
        // visited nodes with nil, to differenciate from nodes in the MDS.
        assert2(bmds.size() == visited.size(), "Wrong bmds input bmds size");
        std::copy(bmds.begin(), bmds.end(), visited.begin());
        std::queue<std::pair<mer_t, mer_t>> queue; // Elements are mer, path length
        for(auto fm : fms) {
            for(mer_t b = 0; b < mer_op_type::alpha; ++b) {
                const auto nm = mer_op_type::nmer(fm, b);
                queue.emplace(nm, 1);
                visited[nm] = nil;
            }
        }

        while(!queue.empty()) {
            const auto elt = queue.front();
            queue.pop();

            unsigned int nb_rc = 0; // Number of right companion of elt in MDS
            for(mer_t b = 0; b < mer_op_type::alpha; ++b) {
                const auto nm = mer_op_type::nmer(elt.first, b);
                if(visited[nm] == yes) {
                    ++nb_rc;
                    continue;
                } else if(visited[nm] == nil) {
                    continue;
                }
                queue.emplace(nm, elt.second + 1);
                visited[nm] = nil;
            }
            if(nb_rc == mer_op_type::alpha) // End of path, it is the shortest
                return elt.second;
        }
        // Should not get there. Didn't find an RF-move, not an MDS.
        assert2(false, "No RF-move found. Not an MDS");
        return 0;
    }
};
#endif // LONGEST_PATH_H_
