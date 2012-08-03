/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshGenerator.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <memory>
#include "BaseLib/BidirectionalMap.h"

#include "MeshLib/Core/Element.h"
#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Core/StructuredMesh.h"

namespace MeshLib
{

class MeshGenerator
{
public:
    static MeshLib::UnstructuredMesh* generateLineMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z);
    static MeshLib::UnstructuredMesh* generateRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z);
    static StructuredMesh<ElementShape::QUAD>* generateStructuredRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z);
    static void generateSubMesh(const MeshLib::IMesh &src, const std::vector<size_t> &list_e, MeshLib::IMesh* &dest, BaseLib::BidirectionalMap<size_t, size_t> &map_global2local);
};

}
