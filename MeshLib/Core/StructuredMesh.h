
#pragma once

#include "IMesh.h"

namespace MeshLib
{

class StructuredMesh : public IMesh
{
private:
    double _origin[3];
    double _length[3];
    double  _unit_length[3];
    size_t  _number_of_nodes_per_dimension[3];
public:
};

}
