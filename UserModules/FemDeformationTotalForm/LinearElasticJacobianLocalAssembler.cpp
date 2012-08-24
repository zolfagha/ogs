/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticJacobianLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "LinearElasticJacobianLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "FemLinearElasticTools.h"
#include "Ogs6FemData.h"


void FemLinearElasticJacobianLocalAssembler::assembly
    (   const NumLib::TimeStep &/*time*/,
        const MeshLib::IElement &e,
        const DiscreteLib::DofEquationIdTable &localDofManager, 
        const NumLib::LocalVector &/*u1*/,
        const NumLib::LocalVector &/*u0*/,
        NumLib::LocalMatrix &localJ )
{
    FemLib::IFiniteElement* fe = _feObjects.getFeObject(e);
    size_t mat_id = e.getGroupID();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    //MaterialLib::PorousMedia* pm = femData->list_pm[mat_id];
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];

    const size_t dim = e.getDimension();
    const size_t n_strain_components = (dim==2 ? 4 : 6);
    const size_t nnodes = e.getNumberOfNodes();
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // set D
    Eigen::MatrixXd matD = Eigen::MatrixXd::Zero(n_strain_components, n_strain_components);
    NumLib::LocalMatrix nv(1,1);
    NumLib::LocalMatrix E(1,1);
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);
    double Lambda, G, K;
    MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
    MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);

    //
    NumLib::LocalMatrix matB(n_strain_components, nnodes*dim);
    NumLib::LocalMatrix matN(dim, nnodes*dim);
    FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
        q->getSamplingPoint(j, gp_x);
        fe->computeBasisFunctions(gp_x);
        fe->getRealCoordinates(real_x);
        double fac = fe->getDetJ() * q->getWeight(j);

        // set N,B
        NumLib::LocalMatrix &N = *fe->getBasisFunction();
        NumLib::LocalMatrix &dN = *fe->getGradBasisFunction();
        setNu_Matrix_byPoint(dim, nnodes, N, matN);
        setB_Matrix_byPoint(dim, nnodes, dN, matB);

        // K += B^T * D * B
        localJ.noalias() += fac * matB.transpose() * matD * matB;
    }
}
