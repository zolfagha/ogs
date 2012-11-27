/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFeObjectContainer.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "LagrangeFeObjectContainer.h"

#include "FemLib/Core/Element/FemElementFactory.h"

namespace FemLib
{

IFiniteElement* LagrangianFeObjectContainer::getFeObject(const MeshLib::IElement &e)
{
    FiniteElementType::type fe_type = getFeType(e.getShapeType(), _order);
    IFiniteElement* fe = FeObjectCachePerFeType::getFeObject(fe_type);
    e.setCurrentOrder(_order);
    fe->configure(*const_cast<MeshLib::IElement*>(&e));
    return fe;
}

IFiniteElement* LagrangianFeObjectContainer::getFeObject(const MeshLib::IElement &e, size_t order)
{
    setPolynomialOrder(order);
    return getFeObject(e);
}


FiniteElementType::type LagrangianFeObjectContainer::getFeType(MeshLib::ElementShape::type ele_type, size_t order)
{
    switch (ele_type)
    {
        case MeshLib::ElementShape::LINE:
            return (order==1) ? FiniteElementType::LINE2 : FiniteElementType::LINE3;
        case MeshLib::ElementShape::QUAD:
            return (order==1) ? FiniteElementType::QUAD4 : FiniteElementType::QUAD9;
        case MeshLib::ElementShape::TRIANGLE:
            return (order==1) ? FiniteElementType::TRI3 : FiniteElementType::TRI6;
        case MeshLib::ElementShape::TETRAHEDRON:
            return (order==1) ? FiniteElementType::TET4 : FiniteElementType::TET10;
        default:
            return FiniteElementType::INVALID;
    }
};

} //end
