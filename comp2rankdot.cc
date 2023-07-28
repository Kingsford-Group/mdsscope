#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <memory>

#include "mer_op.hpp"
#include "mds_op.hpp"
#include "longest_path.hpp"
#include "misc.hpp"
#include "comp2rankdot.hpp"

struct mds_info {
    static size_t total; // Total number of MDSs found. Used to set index

};
size_t mds_info::total = 0;

/* bool operator<(const mds_info& m1, const mds_info& m2) { */
/*     return m1.bmds < m2.bmds; */
/* } */

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef mds_op_type<mer_ops> mds_ops;
typedef longest_path_type<mer_ops> longest_path;

struct fms_index {
    std::vector<mer_t> fms;
    size_t index;
};
typedef std::map<std::vector<tristate_t>, fms_index> layer_type;

std::set<std::vector<tristate_t>> first_layer(const std::vector<tristate_t>& first_bmds) {
    std::set<std::vector<tristate_t>> layer1, layer2;
    layer1.insert(first_bmds);
    std::vector<tristate_t> nbmds;

    // Go back and forth with F-moves and RF-moves, starting from the first MDS,
    // to generate the first 2 full layers. Stop when no new MDS in first layer
    // are discovered.
    bool done = false;
    while(!done) {
        // Do all F-moves from l1 -> l2, then all RF-moves from l2 -> l1. Done
        // if no new MDS added to l1.
        for(const auto& bmds : layer1) {
            for(mer_type fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
                if(mds_ops::has_fm(bmds, fm)) {
                    nbmds = bmds;
                    mds_ops::do_fmove(fm, nbmds);
                    layer2.insert(nbmds);
                }
            }
        }

        done = true;
        for(const auto& bmds : layer2) {
            for(mer_type rfm = 0; rfm < mer_ops::nb_fmoves; ++rfm) {
                if(mds_ops::has_rfm(bmds, rfm)) {
                    nbmds = bmds;
                    mds_ops::do_rfmove(rfm, nbmds);
                    auto it = layer1.insert(nbmds);
                    if(it.second)
                        done = false; // New element in layer1, not done
                }
            }
        }
    }
    return layer1;
}

int main(int argc, char* argv[]) {
    comp2rankdot args(argc, argv);

    std::vector<tristate_t> first_bmds;
    std::vector<mer_t> first_fms, mds;
    mds = mds_from_arg<mer_t>(args.mds_arg);
    mds_ops::from_mds_fms(mds, first_bmds, first_fms);
    auto lp = std::unique_ptr<longest_path>(args.longest_flag ? new longest_path : nullptr);

    layer_type olayer; // Original layer
    {
        const auto mdss = first_layer(first_bmds);
        for(auto& bmds : mdss) {
            std::vector<mer_t> fms;
            for(mer_t fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
                if(mds_ops::has_fm(bmds, fm))
                    fms.push_back(fm);
            }
            fms_index val{std::move(fms), mds_info::total++};
            olayer.emplace(std::move(bmds), std::move(val));
        }
    }

    std::cout << "digraph {\n"
              << "node [shape=circle, style=filled, height=0.2, fixedsize=true];\n"
              << "{ rank = same; \n";
    for(const auto& mdsi : olayer) {
        std::cout << "  n" << mdsi.second.index << " [label=\"\",tooltip=\"" << mdsi.first;
        if(lp)
            std::cout << ':' << lp->longest_path(mdsi.first, mdsi.second.fms);
        std::cout  << "\"];\n";
    }
    std::cout << "}\n";

    // Loop over all layers. stop before the last one so as not to display the
    // first layer multiple times. Save original layer olayer to check
    // connectivity between the last layer and olayer.
    layer_type layer1 = olayer, layer2;
    layer_type* l1 = &layer1;
    layer_type* l2 = &layer2;
    std::vector<tristate_t> nbmds;
    std::vector<mer_t> nfms;
    struct edge_type { size_t n1, n2; mer_t fm; };
    std::vector<edge_type> edges;
    for(size_t fmi = 0; fmi < mer_ops::nb_fmoves - 1; ++fmi, std::swap(l1, l2)) {
        // Find next layer into *l2 and edges between *l1 and *l2
        l2->clear();
        edges.clear();
        for(const auto& mdsi : *l1) {
            // std::cout << "Process " << mdsi.second.index << ": " << mdsi.first << " | " << mdsi.second.fms << '\n';
            const auto& fms = mdsi.second.fms;
            for(auto it = fms.begin(); it != fms.end(); ++it) {
                nbmds = mdsi.first;
                nfms.clear();
                nfms.insert(nfms.end(), fms.begin(), it);
                nfms.insert(nfms.end(), it + 1, fms.end());

                mds_ops::do_fmove(*it, nbmds);
                for(mer_type b = 0; b < mer_ops::alpha; ++b) {
                    const auto nfm = mer_ops::fmove(mer_ops::nmer(*it, b));
                    if(mds_ops::has_fm(nbmds, nfm))
                        nfms.push_back(nfm);
                }
                fms_index val{nfms, 0};
                auto iit  = l2->emplace(nbmds, std::move(val));
                // std::cout << "emplace " << nbmds << " | " << nfms << " -> " << iit.second << '\n';
                if(iit.second) iit.first->second.index = mds_info::total++;
                edges.push_back({mdsi.second.index, iit.first->second.index, *it});
            }
        }

        assert2(!l2->empty(), "Empty layer2");

        // Print the next layer
        std::cout << "{ rank=same;\n";
        for(const auto& mdsi : *l2) {
            std::cout << "  n" << mdsi.second.index << " [label=\"\",tooltip=\"" << mdsi.first;
            if(lp)
                std::cout << ':' << lp->longest_path(mdsi.first, mdsi.second.fms);
            std::cout << "\"];\n";
        }
        std::cout << "}\n";

        // Print edges between the layers
        for(const auto edge : edges) {
            std::cout << "  n" << edge.n1 << " -> " << 'n' << edge.n2 << " [tooltip=\"" << edge.fm << "\"];\n";
        }
    }

    // Every edge from *l1 should now aim to an MDS in olayer. Check that and
    // add edges.
    edges.clear();
    l2 = &olayer;
    for(const auto& mdsi : *l1) {
        for(const auto fm : mdsi.second.fms) {
            nbmds = mdsi.first;
            mds_ops::do_fmove(fm, nbmds);
            const auto it = l2->find(nbmds); // No index. Not keeping this layer
            assert2(it != l2->end(), "Not going to olayer: " << mdsi.first << ' ' << fm << " -> " << nbmds);
            edges.push_back({mdsi.second.index, it->second.index, fm});
        }
    }
    // Print edges between the layers
    for(const auto edge : edges) {
        std::cout << "  n" << edge.n1 << " -> " << 'n' << edge.n2 << " [tooltip=\"" << edge.fm << "\"];\n";
    }

    // End of graph!
    std::cout << "}\n";

    return EXIT_SUCCESS;
}
