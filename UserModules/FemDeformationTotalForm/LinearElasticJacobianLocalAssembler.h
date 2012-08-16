/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticJacobianLocalAssembler.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#pragma once

#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"

class FemLinearElasticJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    explicit FemLinearElasticJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects)
    : _feObjects(feObjects)
    {
    };

    void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ);

private:
    FemLib::LagrangianFeObjectContainer* _feObjects;
};
