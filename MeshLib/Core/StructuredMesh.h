
#pragma once

#include "IMesh.h"

namespace MeshLib
{

class StructuredMesh : public IMesh
{
private:
    GeoLib::Point _origin[3];
    GeoLib::Point _length[3];
    GeoLib::Point  _unit_length[3];
    size_t  _number_of_nodes_per_dimension[3];
public:
};

}
