#ifndef MISC_H_
#define MISC_H_

#include <vector>
#include <unordered_set>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "mer_op.hpp"

template<typename T>
static void append(std::vector<T>& v, const T& x) { v.push_back(x); }
template<typename T>
static void append(std::unordered_set<T>& s, const T& x) { s.insert(x); }

template<typename C>
void mds_from_str(const char* str, C& out) {
    const char* endptr = str;
    while(*endptr) {
        append(out, static_cast<typename C::value_type>(std::strtoul(endptr, (char**)&endptr, 0)));
        endptr += strcspn(endptr, "0123456789");
    }
}

template<typename C>
void mds_parse_to(const std::vector<const char*>& args, C& out) {
    for(auto str : args)
        mds_from_str(str, out);
}

template<typename mer_type>
std::vector<mer_type> mds_from_arg(const std::vector<const char*>& args, bool sort = true) {
    std::vector<mer_type> res;
    mds_parse_to(args, res);
    if(sort)
        std::sort(res.begin(), res.end());
    return res;
}

template<typename C>
void mds_read_to(const char* path_mds, C& out) {
    std::ifstream is(path_mds);
    if(!is.good()) {
		std::ostringstream err;
		err << "Failed to open " << path_mds;
		throw std::runtime_error(err.str());
    }
	uint64_t x = 0;
	while(true) {
		while(is && !std::isdigit(is.peek())) is.get();
		if(!is.good()) break;
		is >> x;
		append(out, static_cast<typename C::value_type>(x));
	}
	if(!is.eof()) [[unlikely]] {
		std::ostringstream err;
		err << "Error while reading " << path_mds;
		throw std::runtime_error(err.str());
	}
}

template<typename mer_type>
std::vector<mer_type> mds_from_file(const char* path_mds, bool sort = true) {
    std::vector<mer_type> res;
    mds_read_to(path_mds, res);
    if(sort)
        std::sort(res.begin(), res.end());
    return res;
}

template<typename C>
void get_mds(const char* path_mds, std::vector<const char*> args, C& res) {
	if(path_mds != nullptr && path_mds[0] != '\0')
		mds_read_to(path_mds, res);
	if(args.size() > 0)
		mds_parse_to(args, res);
}

template<typename C>
inline C get_mds(const char* path_mds, std::vector<const char*> args) {
	C res;
	get_mds(path_mds, args, res);
	return res;
}

#endif // MISC_H_
