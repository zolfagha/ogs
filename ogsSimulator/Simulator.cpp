
#include "Simulator.h"

#include <iostream>

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"
#ifdef USE_LIS
#include "lis.h"
#endif

#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoIO/Ogs4GeoIO.h"
#include "MeshIO/Ogs5MeshIO.h"
#include "FemIO/ogs5/Ogs5FemIO.h"

//#include "ProcessBuilder.h"
#include "FormatterCustom.h"
#include "SimulationInfo.h"
//#include "Ogs5ToOgs6.h"
#include "SimulationProperties.h"

//#include "GeoProcessBuilder.h"

namespace ogs6
{

////////////////////////////////////////////////////////////////////////////////
// Variables
////////////////////////////////////////////////////////////////////////////////
//static FormatterCustom *custom_format;
//static logog::Cout *logogCout;
//static logog::LogFile *logog_file;
//
//static bool isOgsInitCalled = false;
//static bool isOgsExitCalled = false;

//typedef GeoProcessBuilder ProcessBuilder;
////////////////////////////////////////////////////////////////////////////////


void ogsInit(int argc, char* argv[])
{
	if (isOgsInitCalled) return;
	isOgsInitCalled = true;

    LOGOG_INITIALIZE();
    custom_format = new FormatterCustom();
    logogCout = new logog::Cout();
    logogCout->SetFormatter(*custom_format);
    logog_file = NULL;

#ifdef USE_LIS
    lis_initialize(&argc, &argv);
#endif
}

void ogsExit()
{
	if (isOgsExitCalled) return;
	isOgsExitCalled = true;

#ifdef USE_LIS
    lis_finalize();
#endif

	LOGOG_COUT << "exit ogs6." << std::endl;
    BaseLib::releaseObject(custom_format, logogCout, logog_file);
    LOGOG_SHUTDOWN();
}






} //end ogs6
