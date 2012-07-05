
#pragma once
#include <iostream>
#include <fstream>
#include <string>

#include "MeshLib/Core/IMesh.h"

namespace Ogs5MeshIO
{
void readMesh(std::string const& fileName, std::vector<MeshLib::IMesh*> &vec_mesh);
void writeMesh(MeshLib::IMesh const * mesh, std::string const& fileName);
const std::string convertElementType2String(const MeshLib::ElementShape::type t);
MeshLib::ElementShape::type convertString2ElementType(const std::string &s);
}

