
#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include "MeshLib/Core/IMesh.h"

typedef std::pair<std::string, double*> NodalScalarValue;

namespace OutputIOLegacyVtk
{
void WriteAsciiFile(std::string str_file_name, const MeshLib::IMesh &mesh, int time_step_number, double simulation_time, std::vector<NodalScalarValue> &nodalScalarValues);

void WriteVTKHeader(std::fstream &vtk_file, int time_step_number, double simulation_time, std::vector<NodalScalarValue> &nodalScalarValues);

void WriteVTKDataArrays(std::fstream &vtk_file, const MeshLib::IMesh &mesh, std::vector<NodalScalarValue> &nodalScalarValues);
};


