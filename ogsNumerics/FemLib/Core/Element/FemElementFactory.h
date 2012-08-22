/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElementFactory.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <map>

#include "FemLib/Core/Element/FiniteElementType.h"
#include "FemLib/Core/Element/C0IsoparametricElements.h"
#include "FemLib/Core/Element/TRI3CONST.h"

namespace FemLib
{

class FemElementFactory
{
public:
    static IFiniteElement* create(FiniteElementType::type fe_type, MeshLib::IMesh &msh)
    {
        switch (fe_type)
        {
            case FiniteElementType::LINE2:
                return new LINE2(msh);
            case FiniteElementType::LINE3:
                return new LINE3(msh);
            case FiniteElementType::QUAD4:
                return new QUAD4(msh);
            case FiniteElementType::QUAD8:
                return new QUAD8(msh);
            case FiniteElementType::QUAD9:
                return new QUAD9(msh);
            case FiniteElementType::TRI3CONST:
                return new TRI3CONST(msh);
            case FiniteElementType::TRI3:
                return new TRI3(msh);
            case FiniteElementType::TRI6:
                return new TRI6(msh);
            default:
                return 0;
        }
        assert(false);
        return 0;
    }
};

}
