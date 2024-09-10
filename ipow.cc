#include <iostream>
#include <cstdint>

// Print 128 bit long integers. Not very fast. Ignore formatting
std::ostream& operator<<(std::ostream& os, __uint128_t x) {
	static constexpr int buflen = 40;
	char buf[buflen];
	char* ptr = &buf[buflen - 1];
	*ptr = '\0';

    do {
		--ptr;
        *ptr = ((char)(x % 10 + '0'));
        x /= 10;
    } while(x > 0);
	return os << ptr;
}


template<typename T>
// constexpr typename std::enable_if<std::is_integral<T>::value, T>::type
T
ipow(T base, unsigned int exp) {
    T result = 1;
    for(;;) {
        std::cout << exp << ' ' << base << ' ' << result << '\n';
        if(exp & 1) result *= base;
        exp >>= 1;
        if(!exp) break;
        base *= base;
    }

    return result;
}

int main(int argc, char* argv[]) {
	uint64_t a64 = 4;
	__uint128_t a128 = 4;
	int exp = 63;
	std::cout << a64 << "**" << exp << '\n'
			  << ipow(a64, exp) << '\n'
			  << a128 << "**" << exp << '\n'
			  << ipow(a128, exp) << std::endl;

	return 0;
}
