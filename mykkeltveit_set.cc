#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "argparse.hpp"
#include "mykkeltveit.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
#include "pcr_info.hpp"
#include "common.hpp"
#include "dbg.hpp"
#include "longest_path.hpp"

struct MykkeltveiSetArgs : argparse::Args {
    std::optional<uint32_t>& offset_arg = kwarg("o,offset", "Offset in root of unity power").set_default(1);
    bool& all_flag = flag("a,all", "Output all set (both offsets, all combination of repeat PCR)");
    bool& range_flag = flag("r,range", "Find range (smallest/largest longest path)");
    bool& longest_path_flag = flag("l,longest-path", "Compute longest path");
    bool& brute_flag = flag("brute", "Use brute force to create set");

    void welcome() override {
        std::cout <<
            "Create Mykkeltveit\'s set\n\n"
            "Almost identical to Mykkeltveit, but not quite"
            << std::endl;
    }
};

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef pcr_info_type<mer_ops> pcr_info_t;

typedef root_unity_type<mer_ops> root_unity_t;

// Create the Mykkeltveit set. Not optimize at all
std::vector<mer_t> find_mds(const pcr_info_t& pcr_info, unsigned int offset, const root_unity_t& roots) {
    typedef root_unity_t::complex complex;
    const double epsilon = roots.epsilon;

    complex pos, prev_pos, in_pos;
    mer_t prev_m, in_set;
    std::vector<mer_t> mds;
    for(const auto& pcr : pcr_info.pcrs) {
        prev_m = mer_ops::nb_mers;
        prev_pos = complex();
        in_set = mer_ops::nb_mers;
        in_pos = complex();
        for(const auto m : pcr) {
            pos = roots.embed_mer(m, offset);
            // std::cout << "pcr " << (size_t)pcr[0] << " m " << (size_t)m << ' ' << pos <<  " prev " << (size_t)prev_m << ' ' << prev_pos << " in " << (size_t)in_set << ' ' << in_pos << '\n';
            // if(pcr.size() < mer_ops::k || (pos.real() < -epsilon && abs(pos.imag()) < epsilon)) {
            if(abs(pos) < epsilon || (pos.real() < -epsilon && abs(pos.imag()) < epsilon)) {
                in_set = m;
                in_pos = pos;
                break;
            } else if(prev_m != mer_ops::nb_mers && pos.imag() > epsilon && prev_pos.imag() < -epsilon) {
                in_set = prev_m;
                in_pos = prev_pos;
                break;
            }
            prev_pos = pos;
            prev_m = m;
        }
        // Test wrap around PCR if needed
        if(in_set == mer_ops::nb_mers && pos.imag() < -epsilon) {
            // The first k-mer position must be positive, otherwise it is an error
            assert2(roots.embed_mer(pcr[0], offset).imag() > epsilon, "First k-mer should be above x-axis");
            in_set = pcr.back();
            in_pos = pos;
        }
        // std::cout << "end " << (size_t)pcr[0] << " prev " << prev_pos << ' ' << (size_t)prev_m << " in " << (size_t)in_set << ' ' << in_pos << '\n';
        assert2(in_set != mer_ops::nb_mers, "No in set mer for " << mer_ops::alpha << ' ' << mer_ops::k << " PCR " << (size_t)pcr[0]);
        mds.push_back(in_set);
    }
    return mds;
}

std::vector<mer_t> brute_force_mds(const pcr_info_t& pcr_info, unsigned int offset, const root_unity_t& roots) {
    std::vector<mer_t> res;
    for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
        if(roots.in_mykkeltveit_set(m, offset))
            res.push_back(m);
    }
    return res;
}

// Repeat PCR index and offset
struct repeat {
    const std::vector<mer_t>& pcr;
    mer_t index, offset;
};

int main(int argc, char* argv[]) {
    const auto args = argparse::parse<MykkeltveiSetArgs>(argc, argv);

    const pcr_info_t pcr_info;
    root_unity_type<mer_ops> root_unity;

    mer_t offset = *args.offset_arg;
    auto mds = !args.brute_flag ? find_mds(pcr_info, offset, root_unity) : brute_force_mds(pcr_info, offset, root_unity);

    std::vector<repeat> repeat_pcrs;
    if(args.all_flag || args.range_flag) {
        for(mer_t index = 0; index < pcr_info.pcrs.size(); ++index) {
            const auto& pcr = pcr_info.pcrs[index];
            if(pcr.size() > 1 && pcr.size() < mer_ops::k)
                repeat_pcrs.push_back({pcr, index, 0});
        }
    }

    std::vector<mer_t> mds_sorted, min_mds, max_mds;
    size_t min_lp = std::numeric_limits<size_t>::max(), max_lp = 0;
    bool done = false;
    while(!done) {
        if(!args.range_flag) {
            mds_sorted = mds;
            std::sort(mds_sorted.begin(), mds_sorted.end());
            std::cout << joinT<size_t>(mds_sorted, ',');
        }

        if(args.longest_path_flag || args.range_flag) {
            mds_op_type<mer_ops> mds_ops;
            mds_ops.mds2fmoves(mds);
            if(std::find(mds_ops.done.begin(), mds_ops.done.end(), false) != mds_ops.done.end()) {
                if(!args.range_flag)
                    std::cout << " -1"; // Not an MDS
            } else { // Is an MDS
                longest_path_type<mer_ops> lp;
                size_t mds_lp = lp.longest_path(mds);
                if(!args.range_flag) {
                    std::cout << '\t' << mds_lp;
                } else{
                    if(mds_lp < min_lp) {
                        min_mds = mds;
                        min_lp = mds_lp;
                    }
                    if(mds_lp > max_lp) {
                        max_mds = mds;
                        max_lp = mds_lp;
                    }
                }
            }
        }
        if(!args.range_flag)
            std::cout << '\n';

        // if(std::all_of(repeat_pcrs.begin(), repeat_pcrs.end(), [&](const repeat& r) { return r.offset == r.pcr.size(); })) {
        //     done = true;
        //     continue;
        // }
        ssize_t i = repeat_pcrs.size() - 1;
        for( ; i >= 0; --i) {
            auto& r = repeat_pcrs[i];
            ++r.offset;
            if(r.offset < r.pcr.size()) {
                mds[r.index] = r.pcr[r.offset];
                break;
            }
            r.offset = 0;
            mds[r.index] = r.pcr[r.offset];
        }
        if(i < 0) {
            if(!args.offset_arg && args.all_flag && offset == 1) {
                offset = mer_ops::k / 2 + 1;
                mds = find_mds(pcr_info, offset, root_unity);
            } else {
                done = true;
            }
        }
    }

    if(args.range_flag) {
        std::sort(min_mds.begin(), min_mds.end());
        std::sort(max_mds.begin(), max_mds.end());
        std::cout << joinT<size_t>(min_mds, ',') << '\t' << min_lp << '\n'
                  << joinT<size_t>(max_mds, ',') << '\t' << max_lp << '\n';
    }

    return EXIT_SUCCESS;
}
