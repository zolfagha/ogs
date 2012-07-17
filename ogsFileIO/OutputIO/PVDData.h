
#pragma once

#include <string>
#include <vector>


struct PVDData
{
	struct VTK_Info
	{
		double timestep;
		std::string vtk_file;
	};

	std::vector<VTK_Info> vec_dataset;
};
