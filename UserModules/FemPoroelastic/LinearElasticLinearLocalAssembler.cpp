/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearElasticLinearLocalAssembler.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "LinearElasticLinearLocalAssembler.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientLinearEQSLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/Fluid.h"

#include "../FemDeformationTotalForm/FemLinearElasticTools.h"
#include "Ogs6FemData.h"


void FemLinearElasticLinearLocalAssembler::assembly(
        const NumLib::TimeStep &time,  const MeshLib::IElement &e,
        const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n,
        NumLib::LocalEquation &eqs
        )
{
    // parameters need to be passed
    DiscreteLib::DofEquationIdTable* localDofManager = NULL; //TODO
    size_t u_order = 2;
    size_t p_order = 1;

    //
    std::vector<size_t> local_pos_u;
    std::vector<size_t> local_pos_p;
    {
        std::vector<size_t> list_nodeid_u;
        e.getNodeIDList(u_order, list_nodeid_u);
        localDofManager->mapEqsID(0, 0, list_nodeid_u, local_pos_u);
        std::vector<size_t> list_nodeid_p;
        e.getNodeIDList(p_order, list_nodeid_p);
        localDofManager->mapEqsID(0, 0, list_nodeid_p, local_pos_p);

    }
    LocalVectorType u0, p0, u1, p1;
    {
        DiscreteLib::getLocalVector(local_pos_u, local_u_n, u0);
        DiscreteLib::getLocalVector(local_pos_u, local_u_n1, u1);
        DiscreteLib::getLocalVector(local_pos_p, local_u_n, p0);
        DiscreteLib::getLocalVector(local_pos_p, local_u_n1, p1);
    }

    // ------------------------------------------------------------------------
    // Element
    // ------------------------------------------------------------------------
    const size_t dim = e.getDimension();
    const size_t n_strain_components = getNumberOfStrainComponents(dim);
    const size_t nnodes_u = e.getNumberOfNodes(u_order);
    const size_t nnodes_p = e.getNumberOfNodes(p_order);
    const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e.getID());

    // ------------------------------------------------------------------------
    // Transient
    // ------------------------------------------------------------------------
    const double dt = time.getTimeStepSize();

    // ------------------------------------------------------------------------
    // Material (assuming element constant)
    // ------------------------------------------------------------------------
    size_t mat_id = e.getGroupID();
    size_t fluid_id = 0; //TODO
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    MaterialLib::PorousMedia* pm = femData->list_pm[mat_id];
    MaterialLib::Solid *solidphase = femData->list_solid[mat_id];
    MaterialLib::Fluid *fluidphase = femData->list_fluid[fluid_id];

    // solid
    double rho_s = .0;
    solidphase->density->eval(e_pos, rho_s);
    LocalMatrixType De = LocalMatrixType::Zero(n_strain_components, n_strain_components);
    NumLib::LocalMatrix nv(1,1);
    NumLib::LocalMatrix E(1,1);
    solidphase->poisson_ratio->eval(e_pos, nv);
    solidphase->Youngs_modulus->eval(e_pos, E);
    double Lambda, G, K;
    MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
    MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, De);

    // fluid
    double mu = .0;
    fluidphase->dynamic_viscosity->eval(e_pos, mu);
    double rho_f = .0;
    fluidphase->density->eval(e_pos, rho_f);

    // media
    LocalMatrixType k;
    pm->permeability->eval(e_pos, k);
    double n = .0;
    pm->porosity->eval(e_pos, n);
    double s = .0;
    pm->storage->eval(e_pos, s);
    LocalMatrixType k_mu(dim, dim);
    k_mu = k * 1.0/mu;


    // ------------------------------------------------------------------------
    // Body force
    // ------------------------------------------------------------------------
    LocalVectorType body_force(dim);
    body_force *= .0;
    bool hasGravity = false;
    if (hasGravity) {
        body_force[dim-1] = rho_s * 9.81;
    }

    // ------------------------------------------------------------------------
    // Local component assembly
    // ------------------------------------------------------------------------
    LocalMatrixType Kuu = LocalMatrixType::Zero(nnodes_u*dim, nnodes_u*dim);
    LocalMatrixType Cup = LocalMatrixType::Zero(nnodes_u*dim, nnodes_p);
    LocalMatrixType Kpp = LocalMatrixType::Zero(nnodes_p, nnodes_p);
    LocalMatrixType Mpp = LocalMatrixType::Zero(nnodes_p, nnodes_p);
    LocalMatrixType Cpu = LocalMatrixType::Zero(nnodes_p, nnodes_u*dim);
    LocalVectorType Fu = LocalVectorType::Zero(nnodes_u*dim);
    LocalVectorType Fp = LocalVectorType::Zero(nnodes_p);

    // temp matrix
    LocalMatrixType B(n_strain_components, nnodes_u*dim);
    LocalMatrixType Nuvw(dim, nnodes_u*dim);
    const LocalMatrixType m = get_m(dim);

    //
    FemLib::IFiniteElement* fe_u = _feObjects.getFeObject(e, u_order);
    FemLib::IFemNumericalIntegration *q_u = fe_u->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q_u->getNumberOfSamplingPoints(); j++) {
        q_u->getSamplingPoint(j, gp_x);
        fe_u->computeBasisFunctions(gp_x);
        fe_u->getRealCoordinates(real_x);
        double fac_u = fe_u->getDetJ() * q_u->getWeight(j);

        //--- local component ----
        // set N,B
        LocalMatrixType &Nu = *fe_u->getBasisFunction();
        LocalMatrixType &dNu = *fe_u->getGradBasisFunction();
        setNu_Matrix(dim, nnodes_u, Nu, Nuvw);
        setB_Matrix(dim, nnodes_u, dNu, B);
        LocalMatrixType &Np = *fe_u->getBasisFunction();
        LocalMatrixType &dNp = *fe_u->getGradBasisFunction();

        // K_uu += B^T * D * B
        Kuu.noalias() += fac_u * B.transpose() * De * B;

        // C_up += B^T * m * Np
        Cup.noalias() += fac_u * B.transpose() * m * Np;

        // Fu += N^T * b
        Fu.noalias() += fac_u * Nuvw.transpose() * body_force;
    }

    FemLib::IFiniteElement* fe_p = _feObjects.getFeObject(e, p_order);
    FemLib::IFemNumericalIntegration *q_p = fe_p->getIntegrationMethod();
    double gp_x[3], real_x[3];
    for (size_t j=0; j<q_p->getNumberOfSamplingPoints(); j++) {
        q_p->getSamplingPoint(j, gp_x);
        fe_p->computeBasisFunctions(gp_x);
        fe_p->getRealCoordinates(real_x);
        double fac = fe_p->getDetJ() * q_p->getWeight(j);

        //--- local component ----
        // set N,B
        LocalMatrixType &dNu = *fe_p->getGradBasisFunction();
        setB_Matrix(dim, nnodes_u, dNu, B);
        LocalMatrixType &Np = *fe_p->getBasisFunction();
        LocalMatrixType &dNp = *fe_p->getGradBasisFunction();

        // M_pp += Np^T * S * Np
        Mpp.noalias() += fac * Np.transpose() * s * Np;

        // K_pp += dNp^T * K * dNp
        Kpp.noalias() += fac * dNp.transpose() * k_mu * dNp;

        // C_pu += Np^T * m^T * B
        Cpu.noalias() += fac * Np.transpose() * m.transpose() * B;
    }

    // Backward euler
    Fp = 1.0/dt * Mpp * p0 + 1.0/dt * Cpu * u0;

    // ------------------------------------------------------------------------
    // Local equation assembly
    // ------------------------------------------------------------------------
//    LocalMatrixType &localA = *eqs.getA();
//    LocalVectorType &localRHS = *eqs.getRHSAsVec();
    eqs.addAsub(local_pos_u, local_pos_u, Kuu);
    eqs.addAsub(local_pos_u, local_pos_p, Cup, -1.);
    eqs.addAsub(local_pos_p, local_pos_p, Kpp);
    eqs.addAsub(local_pos_p, local_pos_p, Mpp, 1.0/dt);
    eqs.addAsub(local_pos_p, local_pos_u, Cpu, 1.0/dt);

    eqs.addRHSsub(local_pos_u, Fu);
    eqs.addRHSsub(local_pos_p, Fp);
}

