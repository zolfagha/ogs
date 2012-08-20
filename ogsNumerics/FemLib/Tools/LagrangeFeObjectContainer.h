/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFeObjectContainer.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "IFeObjectContainer.h"
#include "FeObjectCachePerFeType.h"

namespace FemLib
{

/**
 * \brief Lagrangian finite element object containers
 */
class LagrangianFeObjectContainer : public FeObjectCachePerFeType, public IFeObjectContainer
{
public:
    explicit LagrangianFeObjectContainer(MeshLib::IMesh &msh) :  FeObjectCachePerFeType(msh), _order(1) {};

    LagrangianFeObjectContainer(const LagrangianFeObjectContainer &src)
    : FeObjectCachePerFeType(*src.getMesh())
    {
        _order = src._order;
    }

    virtual ~LagrangianFeObjectContainer() {};

    void setPolynomialOrder(size_t order) { _order = order; };

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e);

private:
    FiniteElementType::type getFeType(MeshLib::ElementShape::type ele_type, size_t order);

    //DISALLOW_COPY_AND_ASSIGN(LagrangianFeObjectContainer);

private:
    size_t _order;
};



}
