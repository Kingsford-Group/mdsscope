#include <cstdlib>
#include <vector>
#include <set>
#include <map>

#include "mer_op.hpp"
#include "mds_op.hpp"
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
typedef std::map<std::vector<tristate_t>, size_t> layer_type;

layer_type first_layer(const std::vector<tristate_t> first_bmds) {
    layer_type layer1, layer2;
    layer1[first_bmds] = mds_info::total++;
    std::vector<tristate_t> nbmds;

    // Go back and forth with F-moves and RF-moves, starting from the first MDS,
    // to generate the first 2 full layers. Stop when no new MDS in first layer
    // are discovered.
    bool done = false;
    while(!done) {
        // Do all F-moves from l1 -> l2, then all RF-moves from l2 -> l1. Done
        // if no new MDS added to l1.
        for(const auto& mdsi : layer1) {
            for(mer_type fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
                if(mds_ops::has_fm(mdsi.first, fm)) {
                    nbmds = mdsi.first;
                    mds_ops::do_fmove(fm, nbmds);
                    layer2.emplace(nbmds, 0); // No index. Not keeping this layer
                }
            }
        }

        done = true;
        for(const auto& mdsi : layer2) {
            for(mer_type rfm = 0; rfm < mer_ops::nb_fmoves; ++rfm) {
                if(mds_ops::has_rfm(mdsi.first, rfm)) {
                    nbmds = mdsi.first;
                    mds_ops::do_rfmove(rfm, nbmds);
                    auto it = layer1.emplace(nbmds, 0);
                    if(it.second) {
                        it.first->second = mds_info::total++;
                        done = false;
                    }
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

    layer_type olayer = first_layer(first_bmds);

    std::cout << "digraph {\n"
              << "node [shape=circle, height=0.1, fixedsize=true];\n"
              << "{ rank = same; \n";
    for(const auto& mdsi : olayer)
        std::cout << "  n" << mdsi.second << " [label=\"\",tooltip=\"" << mdsi.first << "\"];\n";
    std::cout << "}\n";

    // Loop over all layers. stop before the last one so as not to display the
    // first layer multiple times. Save original layer olayer to check
    // connectivity between the last layer and olayer.
    layer_type layer1 = olayer, layer2;
    layer_type* l1 = &layer1;
    layer_type* l2 = &layer2;
    std::vector<tristate_t> nbmds;
    struct edge_type { size_t n1, n2; mer_t fm; };
    std::vector<edge_type> edges;
    for(size_t fmi = 0; fmi < mer_ops::nb_fmoves - 1; ++fmi, std::swap(l1, l2)) {
        // Find next layer into *l2 and edges between *l1 and *l2
        l2->clear();
        edges.clear();
        for(const auto& mdsi : *l1) {
            for(mer_type fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
                if(mds_ops::has_fm(mdsi.first, fm)) {
                    nbmds = mdsi.first;
                    mds_ops::do_fmove(fm, nbmds);
                    const auto it = l2->emplace(nbmds, 0); // No index. Not keeping this layer
                    if(it.second) it.first->second = mds_info::total++;
                    edges.push_back({mdsi.second, it.first->second, fm});
                }
            }
        }

        // Print the next layer
        std::cout << "{ rank=same;\n";
        for(const auto& mdsi : *l2) {
            std::cout << "  n" << mdsi.second << " [label=\"\",tooltip=\"" << mdsi.first << "\"];\n";
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
        for(mer_type fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
            if(mds_ops::has_fm(mdsi.first, fm)) {
                nbmds = mdsi.first;
                mds_ops::do_fmove(fm, nbmds);
                const auto it = l2->emplace(nbmds, 0); // No index. Not keeping this layer
                if(it.second) {
                    std::cerr << "Not going to olayer: " << mdsi.first << ' ' << fm << " -> " << nbmds << std::endl;
                }
                edges.push_back({mdsi.second, it.first->second, fm});
            }
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
