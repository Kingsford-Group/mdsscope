#include "greedy_mds.hpp"

#include <iostream>
#include <limits>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>
#include "mer_op.hpp"
#include "common.hpp"
#include "pcr_info.hpp"
#include "dbg.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

template<typename mer_op_type>
struct pcr_selection {
    typedef typename mer_op_type::mer_t mer_t;
    const pcr_info_type<mer_op_type>& pcr_info;
    std::vector<mer_t> selection;
    mer_t start_pcr, end_pcr;
    bool done;


    pcr_selection(const pcr_info_type<mer_op_type>& pi)
    : pcr_info(pi)
    , selection(pcr_info.pcrs.size(), 0)
    , start_pcr(0)
    , end_pcr(pcr_info.pcrs.size())
    , done(false)
    {}

    void clear() {
        for(mer_t i = start_pcr; i < end_pcr; ++i)
            selection[i] = 0;
        done = false;
    }

    bool advance() {
        if(done) return false;

        ssize_t i = end_pcr - 1;
        for( ; i >= start_pcr; --i) {
            ++selection[i];
            if(selection[i] < pcr_info.pcrs[i].size()) break;
            selection[i] = 0;
        }
        done = i < start_pcr;

        return true;
    }

    inline bool is_selected(const mer_t mer) const {
        const auto pcr = pcr_info.mer2pcr[mer];
        return mer == pcr_info.pcrs[pcr][selection[pcr]];
    }

    // Copy (sub-)selection from another set.
    void copy(const std::vector<mer_t>& rhs, mer_t from = 0, mer_t len = std::numeric_limits<mer_t>::max()) {
        if(from >= selection.size()) return;
        const auto n = std::min(len, (mer_t)std::min(selection.size() - from, rhs.size()));
        std::copy_n(rhs.begin(), n, selection.begin() + from);
    }
};

template<typename mer_op_type>
struct dfs_dag_type {
    typedef typename mer_op_type::mer_t mer_t;
    std::vector<tristate_t> visiting;
    std::vector<std::pair<mer_t, mer_t>> stack;

    bool is_dag(const pcr_selection<mer_op_type>& selection, mer_t start_pcr = 0, mer_t end_pcr = std::numeric_limits<mer_t>::max()) {
        // Do a DFS on the de Bruijn graph, avoiding the selected nodes, to
        // check if the remaining graph is a DAG
        const pcr_info_type<mer_op_type>& pcr_info = selection.pcr_info;
        visiting.resize(mer_op_type::nb_mers);
        std::fill(visiting.begin(), visiting.end(), nil);
        stack.clear();

        for(mer_t mer = 0; mer < mer_op_type::nb_mers; ++mer) {
            if(visiting[mer] != nil || selection.is_selected(mer) || pcr_info.mer2pcr[mer] < start_pcr || pcr_info.mer2pcr[mer] >= end_pcr)
              continue;
            assert2(stack.empty(), "Stack not empty");

            visiting[mer] = yes;
            stack.emplace_back(mer, 0);

            mer_t m, b;
            while(!stack.empty()) {
                std::tie(m, b) = stack.back();
                if(b >= mer_op_type::alpha) {
                    stack.pop_back();
                    visiting[m] = no;
                    continue;
                }

                ++stack.back().second;
                mer_t nm = mer_op_type::nmer(m, b);
                if(selection.is_selected(nm) || pcr_info.mer2pcr[nm] < start_pcr || pcr_info.mer2pcr[nm] >= end_pcr)
                    continue;
                if(visiting[nm] == yes) return false;
                if(visiting[nm] == nil) {
                    visiting[nm] = yes;
                    stack.emplace_back(nm, 0);
                }
            }
        }
        return true;
    }
};

template<typename mer_op_type>
std::ostream& operator<<(std::ostream& os, const pcr_selection<mer_op_type>& selection) {
    assert2(selection.pcr_info.pcrs.size() == selection.selection.size(), "Selection size differ from # of PCRs");
    if(!selection.pcr_info.pcrs.empty()) {
         os << selection.pcr_info.pcrs[0][selection.selection[0]];
        for(size_t i = 1; i < selection.pcr_info.pcrs.size(); ++i) {
            os << ',' << selection.pcr_info.pcrs[i][selection.selection[i]];
        }
    }
    return os;
}

template<typename mer_op_type>
void thread_worker(pcr_selection<mer_ops>& selection, std::mutex& selection_mtx,
                   mer_t start_pcr,
                   const std::vector<std::vector<typename mer_op_type::mer_t>>& dag_cache,
                   std::mutex& output_mtx) {
    pcr_selection<mer_op_type> nselection(selection.pcr_info);
    dfs_dag_type<mer_op_type> dfs_dag;

    while(true) {
        {
            std::lock_guard<std::mutex> lck(selection_mtx);
            if(!selection.advance()) break;
            nselection.copy(selection.selection, 0, start_pcr);
        }

        if(dfs_dag.is_dag(nselection, 0, start_pcr)) {
            for(const auto& cached : dag_cache) {
                nselection.copy(cached, start_pcr);
                if(dfs_dag.is_dag(nselection)) {
                    std::lock_guard<std::mutex> lck(output_mtx);
                    std::cout << nselection << '\n';
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    greedy_mds args(argc, argv);

    pcr_info_type<mer_ops> pcr_info;
    // std::cerr << "Generated PCR info: " << pcr_info.pcrs.size() << std::endl;
    // for(const auto& pcr : pcr_info.pcrs) {
    //     std::cout << pcr.size() << ':';
    //     for(const auto m : pcr)
    //         std::cout << ' ' << m << '|' << pcr_info.mer2weight[m];
    //     std::cout << '\n';
    // }
    // std::cout << std::flush;
    pcr_selection<mer_ops> selection(pcr_info);

    // Split PCR space in 2: the heavy weights and the lighter ones. Pick just
    // enough PCRs to have no more than some number of possible PCR sets (say <
    // 10^7). Find the subset of those which are acyclic. Only those will be
    // tested later on.
    mer_t start_pcr = pcr_info.pcrs.size();
    size_t nb_pcr_sets = 1;
    while(start_pcr >= 1 && nb_pcr_sets < args.threshold_arg) {
        --start_pcr;
        nb_pcr_sets *= pcr_info.pcrs[start_pcr].size();
    }
    // std::cerr << "Start pcr: " << start_pcr << ' ' << nb_pcr_sets << std::endl;

    selection.start_pcr = start_pcr;
    selection.clear();

    // Fill up cache
    std::vector<std::vector<mer_ops::mer_t>> dag_cache;
    {
        dfs_dag_type<mer_ops> dfs_dag;
        while(true) {
            if(dfs_dag.is_dag(selection, selection.start_pcr))
                dag_cache.emplace_back(selection.selection.cbegin() + start_pcr, selection.selection.cend());
            if(!selection.advance())
                break;
        }
    }

    // std::cerr << "DAG cache: " << dag_cache.size() << '\n';

    // Find all MDSs
    selection.start_pcr = 0;
    selection.end_pcr = start_pcr;
    selection.clear();

    std::vector<std::thread> threads;
    std::mutex selection_mtx, output_mtx;
    const auto nbthreads = args.threads_given ? args.threads_arg : std::thread::hardware_concurrency();
    for(uint64_t i = 0; i < nbthreads; ++i) {
        threads.push_back(std::thread(thread_worker<mer_ops>, std::ref(selection), std::ref(selection_mtx),
                                      start_pcr, std::ref(dag_cache), std::ref(output_mtx)));
    }

    for(auto& th : threads)
        th.join();

    return 0;
}
