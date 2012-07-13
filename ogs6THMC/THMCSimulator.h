
#pragma once

#include <string>
#include <iostream>

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"

#include "BaseLib/CodingTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Options.h"
#include "SimulationInfo.h"
#include "FormatterCustom.h"
#include "GeoProcessBuilder.h"

namespace ogs6
{

/**
 * \brief initialize OGS
 * @param argc
 * @param argv
 */
void ogsInit(int argc, char* argv[]);

/**
 * \brief finalize OGS
 */
void ogsExit();

/**
 * \brief OGS simulator class
 */
class THMCSimulator
{
public:
	typedef GeoProcessBuilder ProcessBuilder;

	THMCSimulator(int argc, char* argv[]);

	~THMCSimulator();

	///
	int execute();

private:
	bool checkInputFiles(const std::string& proj_path);

private:
	SimulationInfo* _sim_info;
};





} //end ogs6
