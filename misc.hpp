#ifndef MISC_H_
#define MISC_H_

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>

#include "mer_op.hpp"

template<typename mer_type>
std::vector<mer_type> mds_from_arg(const std::vector<const char*>& args) {
    std::vector<mer_type> res;
    for(auto str : args) {
        const char* endptr = str;
        while(*endptr) {
            res.push_back(std::strtoul(endptr, (char**)&endptr, 0));
            endptr += strcspn(endptr, "0123456789");
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}


#endif // MISC_H_
