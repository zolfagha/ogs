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

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientResidualLocalAssembler.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"

class FemLinearElasticJacobianLocalAssembler: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    FemLinearElasticJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer* feObjects)
    : _feObjects(feObjects)
    {
    };

    void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ)
    {
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(e);

        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
            q->getSamplingPoint(j, gp_x);
            fe->computeBasisFunctions(gp_x);
            fe->getRealCoordinates(real_x);

//            double k;
//            _pm->hydraulic_conductivity->eval(real_x, k);
//            LocalMatrix matK(1,1);
//
//            //fe->integrateWxN(_pm->storage, localM);
//            fe->integrateDWxDN(j, k, localJ);
        }
    }

private:
    FemLib::LagrangianFeObjectContainer* _feObjects;
};
