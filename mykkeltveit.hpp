#ifndef MYKKELTVEIT_H_
#define MYKKELTVEIT_H_

#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>
#include <iostream>

template<typename mer_ops>
struct root_unity_type {
	typedef std::complex<double> complex;
	typedef typename mer_ops::mer_t mer_t;
	double epsilon = 1e-10;

    std::vector<complex> values;
    root_unity_type() : values(mer_ops::k + 1) {
        for(unsigned int i = 0; i < (mer_ops::k / 2 + 1); ++i) {
            const double theta = 2 * i * M_PI / mer_ops::k;
            values[i] = complex(std::cos(theta), std::sin(theta));
            if(std::abs(values[i].real()) < epsilon)
                values[i].real(0.0);
            if(std::abs(values[i].imag()) < epsilon)
                values[i].imag(0.0);
            if(i > 0 && i < mer_ops::k - i) {
                values[mer_ops::k - i] = std::conj(values[i]);
            }
        }
        values[mer_ops::k] = values[0];
    }

    complex embed_mer(mer_t m, uint32_t offset = 1) const {
        complex pos(0.0, 0.0);
        for(unsigned int i = 0; i < mer_ops::k; ++i, m /= mer_ops::alpha) {
            pos += (complex::value_type)mer_ops::rb(m) * values[(mer_ops::k + offset - 1 - i) % mer_ops::k];
        }
        return pos;
    }

	bool in_mykkeltveit_set(mer_t m, uint32_t offset = 1) const {
		const auto pos = embed_mer(m, offset);
		if((std::abs(pos) < epsilon)) {
			// Choose smallest mer from PCR at the origin.
			for(mer_t nmer = mer_ops::nmer(m); nmer != m; nmer = mer_ops::nmer(nmer))
				if(nmer < m)
					return false;
			return true;
		}

		if(pos.real() < -epsilon && std::abs(pos.imag()) < epsilon) // On the negative x-axis
			 return true;

		// Check if edge to successor cross negative x-axis
		 const mer_t nmer = mer_ops::nmer(m);
		 const auto npos = embed_mer(nmer, offset);
		 return npos.imag() > epsilon && pos.imag() < -epsilon;
	}
};


#endif // MYKKELTVEIT_H_
