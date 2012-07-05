
#pragma once

#include <string>

namespace ogs6
{

void ogsInit(int argc, char* argv[]);

void ogsExit();

class SimulationInfo;

class OgsSimulator
{
public:
	OgsSimulator(int argc, char* argv[]);
	~OgsSimulator();

	int execute();

private:
	bool checkInputFiles(const std::string& proj_path);

private:
	SimulationInfo* _sim_info;

};

} //end ogs6
