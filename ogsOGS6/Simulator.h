
#pragma once

namespace ogs6
{

void ogsInit(int argc, char* argv[]);

void ogsExit();


class Simulator
{
public:
	Simulator(int argc, char* argv[]);
	~Simulator();

	int execute();
};

} //end ogs6
