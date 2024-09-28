#include "argparse.hpp"
#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include <atomic>

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"
#include "mds_op.hpp"
#include "imoves.hpp"
#include "imove_signature.hpp"
#include "file_queue.hpp"
#include "misc.hpp"
#include "backtrace.hpp"

// Only used for debugging.
#ifndef NDEBUG
#include "pcr_info.hpp"
#endif

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;
typedef std::lock_guard<std::mutex> guard_t;

struct TraverseCompArgs : argparse::Args {
    std::string& comps_arg = kwarg("c,comps", "Output file for component");
    std::string& dot_arg = kwarg("d,dot", "Output file for the component graph");
    bool& progress_flag = flag("p,progress", "Display progress");
    uint32_t& threads_arg = kwarg("t,threads", "Thread target (all)").set_default(0);
    std::vector<const char*>& comp_arg = arg("component");

    void welcome() {
        std::cout <<
            "Traverse the component graph\n"
            "Start from one MDS, output all the components in the graph."
            << std::endl;
    }
};

void thread_work(comp_queue<mer_ops>& queue, std::mutex& qlock,
                 signatures_type<mer_ops>& signatures, std::mutex& sig_lock,
                 std::ostream& dot_fd, std::mutex& dot_lock) {
    queue_elt<mer_ops> current, nelt;
    mds_op_type<mer_ops> mds_op;
    imoves_type<mer_ops> imoves_op;
    std::pair<signatures_type<mer_ops>::iterator, bool> ires;
    bool empty = false;
    int retries = 0;

    nelt.fms.resize(mer_ops::nb_mers);

    while(true) {
        {
            guard_t guard(qlock);
            empty = queue.empty();
            if(!empty && !queue.dequeue(current))
                throw std::runtime_error("Failed to dequeue element");
        }
        if(empty) {
            if(retries > 1) break;
            ++retries;
            std::this_thread::sleep_for(std::chrono::seconds(1));
            continue;
        } else {
            retries = 0;
        }
        // std::cout << "current " << current << std::endl;
        mds_op.fromFmoves(current.fms);

        for(const auto im : current.ims) {
            // std::cout << "i-move " << im << std::endl;
            mds_op.traverse_imove(im);
            // assert2(pcr_info.check_bmds(mds_op.nbmds), "Invalid bmds after traversing I-move");
            imoves_op.imoves(mds_op.nbmds, nelt.ims);
            // std::cout << "nfmoves " << mds_op.nfmoves << " imoves " << nelt.ims << std::endl;
            nelt.fms.swap(mds_op.nfmoves);
            {
                guard_t guard(sig_lock);
                ires = signatures.insert(std::make_pair(nelt.ims, signatures.size()));
            }

            // std::cout << "insert " << ires.second << ' ' << ires.first->second << " | " << ires.first->first << std::endl;
            // std::cout << "signatures " << signatures << std::endl;
            nelt.index = ires.first->second;

            // Output dot graph
            // if(current.index < nelt.index) {
            {
                guard_t guard(dot_lock);
                dot_fd << "  n" << current.index
                       << " -> n" << nelt.index
                       << " [label=\"" << im << "\"];\n";
            }
            // }

            if(!ires.second) continue;
            {
                guard_t guard(qlock);
                queue.enqueue(nelt);
            }
        }
    }
}

void progress_thread(const size_t th_target, std::atomic<size_t>& joined,
                     signatures_type<mer_ops>& signatures, std::mutex& sig_lock) {
    size_t size, psize = 0;
    using namespace std::chrono;
    const auto start(system_clock::now());
    while(joined < th_target) {
        {
            guard_t guard(sig_lock);
            size = signatures.size();
        }
        const auto now(system_clock::now());
        const auto diff = size - psize;
        const auto dur = now - start;
        std::cerr << '\r' << size << ' '
                  << diff << ' '
                  << ((double)diff / (dur.count() * system_clock::period::num / system_clock::period::den))
                  << "\033[0k"<< std::flush;
        psize = size;

        std::this_thread::sleep_for(seconds(1));
    }
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    const auto args = argparse::parse<TraverseCompArgs>(argc, argv);

    std::ofstream dot_fd(args.dot_arg);
    if(!dot_fd.good()) {
        std::cerr << "Failed to open " << args.dot_arg << std::endl;
        return EXIT_FAILURE;
    }
    dot_fd << "digraph {\n";
    comp_queue<mer_ops> queue(args.comps_arg.c_str());
    signatures_type<mer_ops> signatures;

    { // Initialize signature set and queue
        imoves_type<mer_ops> imoves_op;
        mds_op_type<mer_ops> mds_op;
        queue_elt<mer_ops> current;
#ifndef NDEBUG
        pcr_info_type<mer_ops> pcr_info;
#endif

        const auto start(mds_from_arg<mer_t>(args.comp_arg));
        assert2(pcr_info.check_mds(start), "Invalid starting MDS");
        current.index = 0;
        current.ims = imoves_op.imoves(start);
        signatures.insert(std::make_pair(current.ims, current.index));

        mds_op.mds2fmoves(start);
        current.fms = mds_op.fmoves;
        queue.enqueue(current);
    }

    std::mutex qlock, sig_lock, dot_lock;
    size_t th_target = args.threads_arg;
    if(th_target == 0) th_target = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    for(size_t i = 0; i < th_target; ++i) {
        threads.emplace_back(thread_work,
                             std::ref(queue), std::ref(qlock),
                             std::ref(signatures), std::ref(sig_lock),
                             std::ref(dot_fd), std::ref(dot_lock));
    }

    std::atomic<size_t> joined(0);
    if(args.progress_flag)
        threads.emplace_back(progress_thread,
                             th_target, std::ref(joined),
                             std::ref(signatures), std::ref(sig_lock));

    for(auto& th : threads) {
        th.join();
        ++joined;
    }

    dot_fd << "}\n";

    return 0;
}
