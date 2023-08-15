#ifndef PCR_INFO_H_
#define PCR_INFO_H_

#include <vector>
#include <limits>
#include <algorithm>

#include "mer_op.hpp"
#include "common.hpp"

template<typename mer_op_type>
struct pcr_info_type {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;
    std::vector<mer_t> mer2weight;
    std::vector<mer_t> mer2pcr;
    std::vector<std::vector<mer_t>> pcrs;

    pcr_info_type()
    : mer2weight(mer_op_t::nb_mers)
    , mer2pcr(mer_op_t::nb_mers, std::numeric_limits<mer_t>::max()) // No info <-> max value
    {
        fill_weights();
        generate_PCRs();
    }

    void fill_weights() {
        for(mer_t m = 0; m < mer_op_type::nb_mers; ++m)
            mer2weight[m] = mer_op_type::weight(m);
    }

    void generate_PCRs() {
        std::vector<mer_t> order(mer_op_type::nb_mers);
        for(mer_t i = 0; i < mer_op_type::nb_mers; ++i)
            order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](mer_t i, mer_t j) -> bool { return mer2weight[i] < mer2weight[j]; });

        for(mer_t i = 0; i < mer_op_type::nb_mers; ++i) {
            const auto start = order[i];
            if(mer2pcr[start] < pcrs.size()) continue; // Already visited

            std::vector<mer_t> new_pcr;
            new_pcr.push_back(start);
            mer2pcr[start] = pcrs.size();
            for(mer_t nmer = mer_op_type::nmer(start); nmer != start; nmer = mer_op_type::nmer(nmer)) {
                new_pcr.push_back(nmer);
                mer2pcr[nmer] = pcrs.size();
            }
            pcrs.emplace_back(new_pcr);
        }
    }

    bool check_mds(const std::vector<mer_t>& mds) {
        std::vector<bool> used_pcrs(pcrs.size(), false);
        if(mds.size() != pcrs.size()) return false;
        for(const auto m : mds) {
            if(m >= mer_op_t::nb_mers) return false;
            const auto pcr = mer2pcr[m];
            if(used_pcrs[pcr]) return false;
            used_pcrs[pcr] = true;
        }
        return true;
    }

    bool check_bmds(const std::vector<tristate_t>& bmds) {
        std::vector<bool> used_pcrs(pcrs.size(), false);
        if(bmds.size() != mer_op_t::nb_mers) return false;
        mer_t count = 0;
        for(mer_t m = 0; m < mer_op_t::nb_mers; ++m) {
            if(bmds[m] != yes) continue;
            ++count;
            const auto pcr = mer2pcr[m];
            if(used_pcrs[pcr]) return false;
            used_pcrs[pcr] = true;
        }
        assert2(count == pcrs.size(), "Wrong number of mers in bmds");
        return true;

    }
};

#endif // PCR_INFO_H_
