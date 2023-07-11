#ifndef FILE_QUEUE_H_
#define FILE_QUEUE_H_

#include <fstream>

#include "mer_op.hpp"
#include "imove_signature.hpp"

template<typename mer_op_type>
struct queue_elt {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;
    typedef imove_sig_type<mer_op_type> imove_sig_t;

    size_t index;
    std::vector<mer_t> fms;
    imove_sig_t ims;
};

template<typename mer_op_type>
std::istream& operator>>(std::istream& is, queue_elt<mer_op_type>& elt) {
    elt.fms.resize(mer_op_type::nb_fmoves);
    elt.ims.clear();

    // Index
    is >> elt.index;
    if(is.peek() != '\t') goto parse_failure;
    if(!is.good()) return is;

    // FMs
    for(size_t i = 0; i < elt.fms.size(); ++i)
        is >> elt.fms[i];
    if(is.peek() != '\t') goto parse_failure;
    if(!is.good()) return is;

    // NB IMs
    size_t nb_ims;
    is >> nb_ims;
    if(is.peek() != '\t') goto parse_failure;
    if(!is.good()) return is;

    // IMs
    elt.ims.resize(nb_ims);
    for(size_t i = 0; i < nb_ims; ++i)
        is >> elt.ims[i];

    return is;

parse_failure:
    is.setstate(std::ios::badbit);
    return is;
}

template<typename mer_op_type>
std::ostream& operator<<(std::ostream& os, const queue_elt<mer_op_type>& elt) {
    // Index
    os << elt.index << '\t';
    if(!os.good()) return os;

    // FMs
    os << elt.fms[0];
    for(size_t i = 1; i < elt.fms.size(); ++i)
        os << ' ' << elt.fms[i];
    if(!os.good()) return os;

    // NB IMs
    os << '\t' << elt.ims.size() << '\t';
    if(!os.good()) return os;

    // IMs
    if(!elt.ims.empty())
        os << elt.ims[0];
    for(size_t i = 1; i < elt.ims.size(); ++i)
        os << ' ' << elt.ims[i];
    return os;
}

// Queue made by reading/writing a file from "both ends".
template<typename mer_op_type>
struct comp_queue {
    typedef mer_op_type mer_op_t;
    typedef typename mer_op_type::mer_t mer_t;
    typedef imove_sig_type<mer_op_type> imove_sig_t;
    typedef queue_elt<mer_op_type> queue_elt_t;

    std::ofstream right; // Must open writing end first
    std::ifstream left;

    comp_queue(const char* path)
    : right(path, std::ios::out | std::ios::trunc)
    , left(path)
    {
        if(!left.good())
            throw std::runtime_error("Failed to open queue file for reading");
        if(!right.good())
            throw std::runtime_error("Failed to open queue file for writing");
    }

    bool empty() { return left.tellg() == right.tellp(); }

    // Read and return first element, advancing by 1 element. Returns false if
    // not successful.
    bool dequeue(queue_elt_t& elt) {
        if(empty()) return false;
        left >> elt;
        if(left.peek() != '\n') return false;
        if(!left.good()) return false;
        left.get(); // Skip new line

        return left.good();
    }

    bool enqueue(const queue_elt_t& elt) {
        right << elt << std::endl; // std::flush;
        return right.good();
    }
};

#endif // FILE_QUEUE_H_