#include <cstdlib>
#include <iostream>

#ifdef __GNUC__
#define G_GNUC_CHECK_VERSION(major, minor) \
    ((__GNUC__ > (major)) || \
     ((__GNUC__ == (major)) && \
      (__GNUC_MINOR__ >= (minor))))
#else
#define G_GNUC_CHECK_VERSION(major, minor) 0
#endif

int main(int argc, char* argv[]) {
	if(argc != 3) {
		std::cerr << "Usage: " << argv[0] << " major minor" << std::endl;
	}
	const int res = G_GNUC_CHECK_VERSION(atoi(argv[1]), atoi(argv[2])) ? EXIT_SUCCESS : EXIT_FAILURE;
	return res;
}
