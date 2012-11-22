/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMmfFemElementFactory.h
 *
 * Created on 2012-11-22 by Norihiro Watanabe
 */

#include "THMmfFemElementFactory.h"

#include <cassert>

#include "FemLib/Core/Element/FiniteElementType.h"
#include "FemLib/Core/Element/FemElementFactory.h"
#include "C0InterfaceElements.h"

namespace THMmf
{

FemLib::IFiniteElement* THMmfFemElementFactory::create(THMmfFiniteElementType::type fe_type, MeshLib::IMesh* msh)
{
    switch (fe_type)
    {
        case THMmfFiniteElementType::LINE2:
        case THMmfFiniteElementType::LINE3:
        case THMmfFiniteElementType::QUAD4:
        case THMmfFiniteElementType::QUAD8:
        case THMmfFiniteElementType::QUAD9:
        case THMmfFiniteElementType::TRI3CONST:
        case THMmfFiniteElementType::TRI3:
        case THMmfFiniteElementType::TRI6:
        case THMmfFiniteElementType::TET4:
        case THMmfFiniteElementType::TET10:
            return FemLib::FemElementFactory::create((FemLib::FiniteElementType::type)fe_type, msh);
        case THMmfFiniteElementType::IE_QUAD4:
            return new IE_QUAD4(msh);
        case THMmfFiniteElementType::IE_TRI3:
            return new IE_TRI3(msh);
        default:
            return 0;
    }
    assert(false);
    return 0;
}

}
