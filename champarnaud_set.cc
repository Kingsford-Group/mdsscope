#include <cstdlib>
#include <vector>
#include <bitset>

#include "divisor.hpp"
#include "misc.hpp"
#include "common.hpp"
#include "dbg.hpp"
#include "champarnaud_set.hpp"

#ifndef K
    #error Must define k-mer length K
#endif

#ifndef ALPHA
    #error Must define alphabet length ALPHA
#endif

#include "mer_op.hpp"

typedef mer_op_type<K, ALPHA> mer_ops;
typedef mer_ops::mer_t mer_t;

template<unsigned alpha_>
struct index_t {
	constexpr static unsigned alpha = alpha_;
	const unsigned m_k;

	std::vector<jflib::divisor64> alphap; // alphap[i] = alpha^i (the divisor)

	index_t(unsigned k_) : m_k(k_) {
		alphap.reserve(m_k+1);
		for(unsigned i = 0; i <= m_k; ++i)
			alphap.emplace_back(ipow(alpha, i));
	}

	inline mer_t base(unsigned k, unsigned i, mer_t m) const {
		assert2(k <= m_k, "Word length at most m_k");
		assert2(i < k, "Base index must be less than k");
		return (m / alphap[k - 1 - i]) % alpha;
	}

	// Get base i (i.e., m[i])
	inline mer_t base(unsigned i, mer_t m) const {
		return base(m_k, i, m);
	}


	// Get subword starting at base i of length l, m[i:l]
	inline mer_t subword(unsigned k, unsigned i, unsigned l, mer_t m) const {
		assert2(k <= m_k, "Word length at most m_k");
		assert2(l > 0, "Length of subword must be non-zero");
		assert2(i < k, "First base of subword must be < k");
		assert2(i + l <= k, "Subword goes past end of word");
		if(i == 0 && l == k) [[unlikely]] return m;
		return (m / alphap[k - i - l]) % alphap[l];
	}

	inline mer_t subword(unsigned i, unsigned l, mer_t m) const {
		return subword(m_k, i, l, m);
	}

    // Standard factorization using Duval's algorithm
	void standard_division(const unsigned k, const mer_t mer, std::vector<unsigned>& res) const {
		unsigned n = 0;
		while(n < k) {
			unsigned i = n + 1, j = n + 2;
			while(true) {
				mer_t bi = base(k, i-1, mer); // Algorithm for 1-based string. Convert to 0-based
				mer_t bj = j <= k ? base(k, j-1, mer) : 0;
				if(j > k || bi > bj) {
					do {
						assert2(j > i, "j must be greater than i or no progress is made");
						n = n + j - i;
						res.push_back(n);
					} while(n < i);
					break;
				} else if(bi < bj) {
					i = n + 1;
					j = j + 1;
				} else { // if(bi = bj)
					i = i + 1;
					j = j + 1;
				}
			}
		}
	}

	inline void standard_division(mer_t mer, std::vector<unsigned>& res) const { return standard_division(m_k, mer, res); }

	inline std::vector<unsigned> standard_division(unsigned k, mer_t mer) const {
		std::vector<unsigned> res;
		standard_division(k, mer, res);
		return res;
	}

	inline std::vector<unsigned> standard_division(mer_t mer) const { return standard_division(m_k, mer); }

    // Lyndon word test. Based on the standard factorization algorithm: the only
    // element returned should be k itself.
	bool is_lyndon(const unsigned k, const mer_t mer) const {
		unsigned n = 0;
		while(n < k) {
			unsigned i = n + 1, j = n + 2;
			while(true) {
				mer_t bi = base(k, i-1, mer); // Algorithm for 1-based string. Convert to 0-based
				mer_t bj = j <= k ? base(k, j-1, mer) : 0;
				if(j > k || bi > bj) {
					assert2(j > i, "j must be greater than i or no progress is made");
					n = n + j - i;
					return n == k;
				} else if(bi < bj) {
					i = n + 1;
					j = j + 1;
				} else { // if(bi = bj)
					i = i + 1;
					j = j + 1;
				}
			}
		}
		return false;
	}

	bool is_lyndon(const mer_t mer) const { return is_lyndon(m_k, mer); }

	// Do the "main division" of mer. I.e., mer = l^n u, where l is the shortest
	// such Lyndon word and |u| < |l|. Returns the length of the Lyndon word |l|
	// and |l|^n.
	std::pair<unsigned, unsigned> parts(const unsigned k, const mer_t mer) const {
		assert2(k <= m_k, "k must be less than m_k");
		for(unsigned i = 1; i < k; ++i) { // Check division by a Lyndon word of length 1
			mer_t l = subword(0, i, mer); // Take subword of k-mer mer
			if(is_lyndon(i, l)) { // Check if it is a Lyndon word of length i
				bool ispower = true; // Check if this Lyndon word is repeated in mer, i.e., mer = l^n u
				unsigned j = i;
				for( ; j + i <= k && ispower; j += i)
					ispower = l == subword(j, i, mer);
				if(ispower)
					return std::make_pair(i, j);
			}
		}
		return std::make_pair(k, k);
	}

	inline std::pair<unsigned, unsigned> parts(const mer_t mer) const { return parts(m_k, mer); }

	// Returns the reversed k-mers with the rest followed by the main part.
	mer_t reversed(const unsigned k, const mer_t mer) const {
		auto div = parts(mer);
		if(div.second == k) return mer;
		return subword(0, div.second, mer) + subword(div.second, k - div.second, mer) * alphap[div.second].d();
	}

	inline mer_t reversed(mer_t mer) const { return reversed(m_k, mer); }

};

template<typename Fn>
void enumerate_champarnaud(Fn fn) {
    index_t<ALPHA> index(K);

	std::vector<bool> visited(mer_ops::nb_mers);
	for(mer_t m = 0; m < mer_ops::nb_mers; ++m) {
		if(visited[m]) continue;
		for(mer_t cur = m; !visited[cur]; cur = mer_ops::nmer(cur)) // Mark PCR as visited
			visited[cur] = true;
		fn(index.reversed(m));
	}

}

int main(int argc, char* argv[]) {
	champarnaud_set args(argc, argv);

    if(args.notsorted_flag) {
        bool first = true;
        auto fn = [&first](mer_t m) {
            if(!first) {
                std::cout << ',';
            } else {
                first = false;
            }
            std::cout << (size_t)m;
        };
        enumerate_champarnaud(fn);
    } else {
        std::vector<mer_t> set;
        auto fn = [&set](mer_t m) { set.push_back(m); };
        enumerate_champarnaud(fn);
        std::sort(set.begin(), set.end());
        std::cout << joinT<size_t>(set, ',');
    }

    std::cout << '\n';

	return EXIT_SUCCESS;
}
