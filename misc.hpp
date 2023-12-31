#ifndef MISC_H_
#define MISC_H_

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>

#include "mer_op.hpp"

template<typename mer_type>
std::vector<mer_type> mds_from_arg(const std::vector<const char*>& args, bool sort = true) {
    std::vector<mer_type> res;
    for(auto str : args) {
        const char* endptr = str;
        while(*endptr) {
            res.push_back(std::strtoul(endptr, (char**)&endptr, 0));
            endptr += strcspn(endptr, "0123456789");
        }
    }
    if(sort)
        std::sort(res.begin(), res.end());
    return res;
}

template<typename mer_t>
std::vector<mer_t> read_mds_from_file(const char* path_mds) {
    std::ifstream is(path_mds);
    if(!is.good()) {
        std::cerr << "Failed to open " << path_mds << std::endl;
        exit(1);
    }

    std::string str(std::istreambuf_iterator<char>{is}, {});
    std::vector<const char*> vs{str.c_str()};
    return mds_from_arg<mer_t>(vs);
}


#endif // MISC_H_
