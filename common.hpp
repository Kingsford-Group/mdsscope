#ifndef COMMON_H_
#define COMMON_H_

#include <cstdint>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>

// State of a mer. Either it is nil (unknown), no (absent), yes (present) or
// blocked (should not take part in an F-move).
enum tristate { nil, no, yes, blocked };
typedef uint8_t tristate_t;

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& mds) {
    for(size_t i = 0; i < mds.size(); ++i) {
        if(i) os << ' ';
        os << mds[i];
    }
    return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const std::vector<tristate_t>& mds);

template<typename C>
std::ostream& print_set(std::ostream& os, const C& set) {
    if(!set.empty()) {
        auto it = set.cbegin();
        os << *it;
        for(++it; it != set.cend(); ++it)
            os << ',' << *it;
    }
    return os;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& set) { return print_set(os, set); }
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::set<T>& set) { return print_set(os, set); }

template<typename T, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T, V>& map) {
    if(!map.empty()) {
        auto it = map.cbegin();
        os << '<' << it->first << '>' << it->second;
        for(++it; it != map.cend(); ++it)
            os << '<' << it->first << '>' << it->second;
    }
    return os;
}

template<typename C>
struct join {
    const C& c;
    const char s;

    join(const C& cont, char sep)
        : c(cont)
        , s(sep)
    {}
};

template<typename C>
std::ostream& operator<<(std::ostream& os, const join<C>& j) {
    auto it = j.c.begin();
    if(it != j.c.end()) {
        os << *it;
        for(++it; it != j.c.end(); ++it)
            os << j.s << *it;
    }
    return os;
}

#endif // COMMON_H_
