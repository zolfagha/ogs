/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMmfFemElementCatalogBuilder.h
 *
 * Created on 2012-11-22 by Norihiro Watanabe
 */

#pragma once

#include "FemLib/Tools/FemElementCatalog.h"
#include "FemLib/Tools/MeshElementGroupAndShapeToFemElementType.h"
#include "FemLib/Core/Element/C0IsoparametricElements.h"
#include "C0InterfaceElements.h"
#include "THMmfFiniteElementType.h"

namespace THMmf
{

class THMmfFemElementCatalogBuilder
{
public:
    static void construct(const std::vector<size_t> &vec_m_group_id, const std::vector<size_t> &vec_ie_group_id, FemLib::FemElementCatalog &feCatalog, FemLib::MeshElementGroupAndShapeToFemElementType &shape2fetype)
    {
        feCatalog.registerFeType<FemLib::LINE2>(THMmfFiniteElementType::LINE2);
        feCatalog.registerFeType<FemLib::LINE3>(THMmfFiniteElementType::LINE3);
        feCatalog.registerFeType<FemLib::QUAD4>(THMmfFiniteElementType::QUAD4);
        feCatalog.registerFeType<FemLib::QUAD8>(THMmfFiniteElementType::QUAD8);
        feCatalog.registerFeType<FemLib::QUAD9>(THMmfFiniteElementType::QUAD9);
        feCatalog.registerFeType<FemLib::TRI3>(THMmfFiniteElementType::TRI3);
        feCatalog.registerFeType<FemLib::TRI6>(THMmfFiniteElementType::TRI6);
        feCatalog.registerFeType<FemLib::TET4>(THMmfFiniteElementType::TET4);
        feCatalog.registerFeType<FemLib::TET10>(THMmfFiniteElementType::TET10);
        feCatalog.registerFeType<IE_TRI3>(THMmfFiniteElementType::IE_TRI3);
        feCatalog.registerFeType<IE_QUAD4>(THMmfFiniteElementType::IE_QUAD4);

        for (size_t i=0; i<vec_m_group_id.size(); ++i) {
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::LINE, 1, THMmfFiniteElementType::LINE2);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::LINE, 2, THMmfFiniteElementType::LINE3);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::QUAD, 1, THMmfFiniteElementType::QUAD4);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::QUAD, 2, THMmfFiniteElementType::QUAD9);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::TRIANGLE, 1, THMmfFiniteElementType::TRI3);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::TRIANGLE, 2, THMmfFiniteElementType::TRI6);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::TETRAHEDRON, 1, THMmfFiniteElementType::TET4);
            shape2fetype.addFeType(vec_m_group_id[i], MeshLib::ElementShape::TETRAHEDRON, 2, THMmfFiniteElementType::TET10);
        }
        for (size_t i=0; i<vec_ie_group_id.size(); ++i) {
            shape2fetype.addFeType(vec_ie_group_id[i], MeshLib::ElementShape::TRIANGLE, 1, THMmfFiniteElementType::IE_TRI3);
            shape2fetype.addFeType(vec_ie_group_id[i], MeshLib::ElementShape::QUAD, 1, THMmfFiniteElementType::IE_QUAD4);
        }
    }
};

}
