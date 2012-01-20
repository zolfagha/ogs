
#pragma once

#include "IMesh.h"

namespace MeshLib
{

template<typename Tpos>
class StructuredMesh : public IMesh<Tpos>
{
private:
    Tpos _origin[3];
    Tpos _length[3];
    Tpos  _unit_length[3];
    size_t  _number_of_nodes_per_dimension[3];
public:
};

}
