
#pragma once

#include <iostream>
#include <fstream>
#include <string>

namespace MeshLib
{
class IMesh;
}

namespace MeshIOLegacyVtk
{
    void WriteAsciiFile(const std::string &str_file_name, const MeshLib::IMesh &mesh);

    void WriteVTKHeader(std::fstream &vtk_file);

    void WriteVTKPoints(std::fstream &vtk_file, const MeshLib::IMesh &mesh);

    void WriteVTKCells(std::fstream &vtk_file, const MeshLib::IMesh &mesh);

    void WriteVTKCellDataHeader(std::fstream &vtk_file, const MeshLib::IMesh &mesh);

    void WriteVTKMaterialID(std::fstream &vtk_file, const MeshLib::IMesh &mesh);
};





