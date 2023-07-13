#include "common.hpp"

template<>
std::ostream& operator<<(std::ostream& os, const std::vector<tristate_t>& mds) {
    bool notfirst = false;
    for(size_t i = 0; i < mds.size(); ++i) {
        if(mds[i] == yes) {
            if(notfirst) {
                os << ' ';
            } else {
                notfirst = true;
            }
            os << i;
        }
    }
    return os;
}
