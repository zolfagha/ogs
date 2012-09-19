/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GWJacobianLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "HeadBasedGWJacobianLocalAssembler.h"

#include "NumLib/Function/TXFunction.h"
#include "Ogs6FemData.h"

void HeadBasedGWJacobianLocalAssembler::assembly(const NumLib::TimeStep &ts, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ)
{
    FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
    size_t mat_id = e.getGroupID();
    MaterialLib::PorousMedia* pm = Ogs6FemData::getInstance()->list_pm[mat_id];
    NumLib::LocalMatrix localM = NumLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());
    NumLib::LocalMatrix localK = NumLib::LocalMatrix::Zero(localJ.rows(), localJ.cols());

    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);

        NumLib::LocalMatrix k;
        pm->hydraulic_conductivity->eval(real_x, k);
        NumLib::LocalMatrix s;
        pm->storage->eval(real_x, s);

        fe->integrateWxN(j, s, localM);
        fe->integrateDWxDN(j, k, localK);
    }

    double euler_theta = 1.0; //TODO where to get this parameter?
    localJ.noalias() = 1./ts.getTimeStepSize() * localM + euler_theta *localK;
}

