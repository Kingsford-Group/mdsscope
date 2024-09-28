#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <fstream>

#include "argparse.hpp"
#include "mer_op.hpp"
#include "mds_op.hpp"
#include "longest_path.hpp"
#include "misc.hpp"
#include "common.hpp"

struct Comp2RankdotArgs : argparse::Args {
    bool& longest_flag = flag("l,longest", "Annotate with longest remaining path");
    std::string& output_arg = kwarg("o,output", "Dot file output").set_default("/dev/stdout");
    bool& progress_flag = flag("p,progress", "Show progress");
    std::vector<const char*>& mds_arg = arg("MDS").set_default("");

    void welcome() override {
        std::cout <<
            "Generate dot file for 1 component, with rank for proper plotting\n\n"
            "From 1 MDS, do a BFS and layout each layer with same rank"
            << std::endl;
    }
};

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

std::set<std::vector<tristate_t>> first_layer(const std::vector<tristate_t>& first_bmds, bool progress) {
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
        if(progress)
            std::cerr << '\r' << layer1.size() << ' ' << layer2.size() << std::flush;
        for(const auto& bmds : layer1) {
            for(mer_t fm = 0; fm < mer_ops::nb_fmoves; ++fm) {
                if(mds_ops::has_fm(bmds, fm)) {
                    nbmds = bmds;
                    mds_ops::do_fmove(fm, nbmds);
                    layer2.insert(nbmds);
                }
            }
        }

        done = true;
        if(progress)
            std::cerr << '\r' << layer1.size() << ' ' << layer2.size() << std::flush;
        for(const auto& bmds : layer2) {
            for(mer_t rfm = 0; rfm < mer_ops::nb_fmoves; ++rfm) {
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

void update_ranges(std::vector<std::pair<mer_t, mer_t>>& fm_ranges, mer_t fmi) {
    if(fm_ranges.empty() || (mer_t)(fm_ranges.back().second + 1) < fmi) {
        fm_ranges.emplace_back(fmi, fmi);
    } else if(fm_ranges.back().second != fmi) {
        fm_ranges.back().second = fmi;
    }
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    const auto args = argparse::parse<Comp2RankdotArgs>(argc, argv);

    std::vector<tristate_t> first_bmds;
    std::vector<mer_t> first_fms, mds;
    mds = mds_from_arg<mer_t>(args.mds_arg);
    mds_ops::from_mds_fms(mds, first_bmds, first_fms);
    auto lp = std::unique_ptr<longest_path>(args.longest_flag ? new longest_path : nullptr);
    std::pair<mer_t, mer_t> lprange{std::numeric_limits<mer_t>::max(), 0}, llprange; // Longest path ranges (global and per layer)
    std::pair<mer_t, mer_t> sprange(std::numeric_limits<mer_t>::max(), 0), lsprange; // Shortest path ranges
    std::pair<mer_t, mer_t> width_range; // Min and max width of a layer
    size_t total_mds = 0; // Total number of MDSs seen


    std::ofstream dot_fd;

    if(args.output_arg.size() > 0) {
        dot_fd.open(args.output_arg);
        if(!dot_fd.good()) {
            std::cerr << "Failed to open dot output file '" << args.output_arg << '\'' << std::endl;
            return EXIT_FAILURE;
        }
    }

    if(args.progress_flag)
        std::cerr << "Initialize" << std::flush;

    layer_type olayer; // Original layer
    {
        const auto mdss = first_layer(first_bmds, args.progress_flag);
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

    width_range.first = width_range.second = olayer.size();
    total_mds = olayer.size();

    if(dot_fd.is_open()) {
        dot_fd << "digraph {\n"
               << "node [shape=circle, style=filled, height=0.2, fixedsize=true];\n"
               << "{ rank = same; \n";
        for(const auto& mdsi : olayer) {
            dot_fd << "  n" << mdsi.second.index << " [label=\"\",tooltip=\"" << mdsi.first;
            if(lp) {
                const auto lpl = lp->longest_path(mdsi.first, mdsi.second.fms);
                dot_fd << ':' << (uint64_t)lpl;
                lprange.first = std::min(lprange.first, lpl);
                lprange.second = std::max(lprange.second, lpl);
            }
            dot_fd  << "\"];\n";
        }
        dot_fd << "}\n";
    }

    // Loop over all layers. stop before the last one so as not to display the
    // first layer multiple times. Save original layer olayer to check
    // connectivity between the last layer and olayer.
    layer_type layer1 = olayer, layer2;
    layer_type* l1 = &layer1;
    layer_type* l2 = &layer2;
    std::vector<std::vector<std::pair<mer_t, mer_t>>> fms_ranges(mer_ops::nb_fmoves); // Ranges where F-move is used
    std::vector<tristate_t> nbmds;
    std::vector<mer_t> nfms;
    struct edge_type { size_t n1, n2; mer_t fm; };
    std::vector<edge_type> edges;
    bool is_mds = true;
    for(size_t fmi = 0; fmi + 1 < mer_ops::nb_fmoves; ++fmi, std::swap(l1, l2)) {
        if(args.progress_flag)
            std::cerr << '\r' << fmi << ' ' << total_mds << ' ' << (size_t)width_range.first << ':' << (size_t)width_range.second << std::flush;
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
                for(mer_t b = 0; b < mer_ops::alpha; ++b) {
                    const auto nfm = mer_ops::fmove(mer_ops::nmer(*it, b));
                    if(mds_ops::has_fm(nbmds, nfm))
                        nfms.push_back(nfm);
                }
                fms_index val{nfms, 0};
                auto iit  = l2->emplace(nbmds, std::move(val));
                if(iit.second) iit.first->second.index = mds_info::total++;
                edges.push_back({mdsi.second.index, iit.first->second.index, *it});
                update_ranges(fms_ranges[*it], fmi);
            }
        }

        if(l2->empty()) {
            is_mds = false;
            break;
        }

        // Print the next layer
        llprange.first = std::numeric_limits<mer_t>::max();
        llprange.second = 0;
        lsprange.first = std::numeric_limits<mer_t>::max();
        lsprange.second = 0;
        if(l2->size() < width_range.first)
            width_range.first = l2->size();
        if(l2->size() > width_range.second)
            width_range.second = l2->size();
        total_mds += l2->size();

        if(dot_fd.is_open()) {
            dot_fd << "{ rank=same;\n";
            for(const auto& mdsi : *l2) {
                dot_fd << "  n" << mdsi.second.index << " [label=\"\",tooltip=\"" << mdsi.first;
                if(lp) {
                    const auto spl = lp->shortest_path(mdsi.first, mdsi.second.fms);
                    const auto lpl = lp->longest_path(mdsi.first, mdsi.second.fms);
                    dot_fd << ':' << (uint64_t)spl << ':' << (uint64_t)lpl;
                    lprange.first = std::min(lprange.first, lpl);
                    lprange.second = std::max(lprange.second, lpl);
                    llprange.first = std::min(llprange.first, lpl);
                    llprange.second = std::max(llprange.second, lpl);
                    sprange.first = std::min(sprange.first, spl);
                    sprange.second = std::max(sprange.second, spl);
                    lsprange.first = std::min(lsprange.first, spl);
                    lsprange.second = std::max(lsprange.second, spl);

                }
                dot_fd << "\"];\n";
            }
            dot_fd << "}\n";

            // Print edges between the layers
            for(const auto edge : edges) {
                dot_fd << "  n" << edge.n1 << " -> " << 'n' << edge.n2 << " [tooltip=\"" << (size_t)edge.fm << "\"];\n";
            }
        }
    }

    // Every edge from *l1 should now aim to an MDS in olayer. Check that and
    // add edges.
    if(is_mds && dot_fd.is_open()) {
        edges.clear();
        l2 = &olayer;
        for(const auto& mdsi : *l1) {
            for(const auto fm : mdsi.second.fms) {
                nbmds = mdsi.first;
                mds_ops::do_fmove(fm, nbmds);
                const auto it = l2->find(nbmds); // No index. Not keeping this layer
                assert2(it != l2->end(), "Not going to olayer: " << mdsi.first << ' ' << fm << " -> " << nbmds);
                edges.push_back({mdsi.second.index, it->second.index, fm});
                update_ranges(fms_ranges[fm], mer_ops::nb_fmoves - 1);
            }
        }
        // Print edges between the layers
        for(const auto edge : edges) {
            dot_fd << "  n" << edge.n1 << " -> " << 'n' << edge.n2 << " [tooltip=\"" << (size_t)edge.fm << "\"];\n";
        }
    }

    // End of graph!
    if(dot_fd.is_open())
        dot_fd << "}\n";

    if(args.progress_flag)
        std::cerr << std::endl;

    if(lp)
        std::cerr << "ranges: " << (size_t)sprange.first << ' ' << (size_t)sprange.second << ' ' << (size_t)lprange.first << ' ' << (size_t)lprange.second << "\n";
    std::cerr << "mds: " << is_mds << '\n';
    std::cerr << "width: " << (size_t)width_range.first << ' ' << (size_t)width_range.second << '\n';
    std::cerr << "total: " << total_mds << '\n';

    // bool split_ranges = false;
    // for(auto it = fms_ranges.cbegin(); it != fms_ranges.cend(); ++it) {
    //     const size_t fm = it - fms_ranges.cbegin();
    //     assert2(it->size() > 0, "F-move not used in graph " << fm);
    //     size_t tr = 0; // sum of ranges
    //     for(const auto& r : *it)
    //         tr += r.second - r.first + 1;
    //     auto nb = it->size();
    //     std::cout << fm << ' ' << nb;
    //     for(const auto& r : *it)
    //         std::cout << ' ' << (size_t)r.first << '-' << (size_t)r.second;
    //     std::cout << '\n';
    //     if(nb > 1 && it->front().first == 0 && it->back().second == mer_ops::nb_fmoves - 1) // Wrap around? One less range really
    //         --nb;
    //     if(nb != 1 || tr == mer_ops::nb_fmoves) {
    //         std::cout << fm << ' ' << nb << ' ' << tr << '\n';
    //         // split_ranges = true;
    //     }
    // }
    // for(mer_t fm1 = 0; fm1 < mer_ops::nb_fmoves; ++fm1) {
    //     std::cout << (size_t)fm1 << ':';
    //     for(mer_t fm2 = 0; fm2 < mer_ops::nb_fmoves; ++fm2) {
    //         if(fm1 == fm2) continue;
    //         bool overlap = false;
    //         for(const auto& r1 : fms_ranges[fm1]) {
    //             for(const auto& r2 : fms_ranges[fm2]) {
    //                 overlap = overlap || (r1.second >= r2.first && r1.first <= r2.second);
    //             }
    //         }
    //         if(!overlap)
    //             std::cout << ' ' << (size_t)fm2;
    //     }
    //     std::cout << '\n';
    // }

    return EXIT_SUCCESS;
}
