
#pragma once

#include <string>
#include <vector>


struct PVDData
{
	struct
	{
		double timestep;
		std::string vtk_file;
	} VTK_Info;

	std::vector<VTK_Info> vec_dataset;
};
