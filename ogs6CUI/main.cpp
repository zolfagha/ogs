/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file main.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include <iostream>
#include <exception>
#include "THMCSimulator.h"

int main ( int argc, char *argv[] )
{
    argc = 5;
    argv = (char**) malloc((argc+1) * sizeof(char *));
    argv[0] = "ogs6";
    argv[1] = "-i";
    argv[2] = "E:\\3.Task\\20120719_ogs6test\\Mass\\q_quad";
    argv[3] = "-o";
    argv[4] = "E:\\3.Task\\20120719_ogs6test\\Mass\\ogs6";

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
