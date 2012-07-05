
#include <exception>
#include "ogs6.h"

int main ( int argc, char *argv[] )
{
	int returncode = 0;
	try {
//		ogs6::ogsInit(argc, argv);

		ogs6::Simulator sim(argc, argv);
		returncode = sim.execute();

	} catch (char* e) {
		std::cerr << e << std::endl;
	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
	} catch (...) {
		std::cerr << "Unknown exception occurred!" << std::endl;
	}

//	ogs6::ogsExit(); // finalize must be called!

	return returncode;
}
