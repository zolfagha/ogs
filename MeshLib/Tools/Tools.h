
#pragma once

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "MeshLib/Core/IMesh.h"


namespace MeshLib
{

template<typename Tpos>
void findNodesOnGeometry(const IMesh<Tpos> *msh, const GeoLib::Polyline *poly, std::vector<Node<Tpos>*> *vec_nodes)
{
  throw std::exception("The method or operation is not implemented.");
};

template<typename Tpos>
void findElementsOnGeometry(const IMesh<Tpos> *msh, const GeoLib::Polyline *poly, std::vector<size_t> *vec_nodes)
{
  throw std::exception("The method or operation is not implemented.");
};

}