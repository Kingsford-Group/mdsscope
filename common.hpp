#ifndef COMMON_H_
#define COMMON_H_

#include <iostream>
#include <vector>
#include <unordered_set>

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

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& set) {
    if(!set.empty()) {
        auto it = set.cbegin();
        os << *it;
        for(++it; it != set.cend(); ++it)
            os << ',' << *it;
    }
    return os;
}

#endif // COMMON_H_
