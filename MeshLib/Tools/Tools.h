
#pragma once

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "MeshLib/Core/IMesh.h"


namespace MeshLib
{

template<typename Tpos>
void findNodesOnGeometry(const IMesh<Tpos> *msh, const GeoLib::Polyline *poly, std::vector<Node<Tpos>*> *vec_nodes)
{
    
};


}