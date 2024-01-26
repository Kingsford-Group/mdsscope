#include "create_seed.hpp"
#include "random_seed.hpp"

int main(int argc, char* argv[]) {
	create_seed args(argc, argv);
	seeded_prg<std::mt19937_64>(args.seed_arg, nullptr);
	return EXIT_SUCCESS;
}
