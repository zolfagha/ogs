
#include <iostream>
#include <exception>
#include "THMCSimulator.h"

int main ( int argc, char *argv[] )
{
    int returncode = 0;
    try {
        ogs6::THMCSimulator sim(argc, argv);
        returncode = sim.execute();
    } catch (char* e) {
        std::cerr << e << std::endl;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception occurred!" << std::endl;
    }

    return returncode;
}
