#include "argparse.hpp"
#include "random_seed.hpp"

struct CreateSeedArgs : argparse::Args {
	std::string& seed_arg = arg("Output seed file");

	void welcome() {
		std::cout << "Create a random seed file" << std::endl;
	}
};

int main(int argc, char* argv[]) {
	const auto args = argparse::parse<CreateSeedArgs>(argc, argv);
	seeded_prg<std::mt19937_64>(args.seed_arg.c_str(), nullptr);
	return EXIT_SUCCESS;
}
