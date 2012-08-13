/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StressStrain.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "StressStrain.h"

#include "logog.hpp"

#include "Ogs6FemData.h"
#include "FemLinearElasticTools.h"


size_t FunctionStressStrain::getNumberOfStrainComponents() const
{
    MeshLib::IMesh* msh = _dis->getMesh();
    const size_t dim = msh->getDimension();
    return (dim==2 ? 4 : 6);
}

bool FunctionStressStrain::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    size_t msh_id = option.getOption<size_t>("MeshID");
    _dis = femData->list_dis_sys[msh_id];
    const size_t n_strain_components = getNumberOfStrainComponents();

    // create strain, stress vectors
    _strain = new FemLib::FEMIntegrationPointFunctionVector(*_dis);
    _stress = new FemLib::FEMIntegrationPointFunctionVector(*_dis);

    // set initial output
    OutputVariableInfo var1("STRAIN", OutputVariableInfo::Element, OutputVariableInfo::Real, n_strain_components, _strain);
    femData->outController.setOutput(var1.name, var1);
    OutputVariableInfo var2("STRESS", OutputVariableInfo::Element, OutputVariableInfo::Real, n_strain_components, _stress);
    femData->outController.setOutput(var2.name, var2);

    // initial output parameter
    setOutput(Strain, _strain);
    setOutput(Stress, _stress);

    return true;
}

void FunctionStressStrain::accept(const NumLib::TimeStep &/*time*/)
{
    //std::cout << "Velocity=" << std::endl;
    //_vel->printout();
    //update data for output
    const size_t n_strain_components = getNumberOfStrainComponents();
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var1("STRAIN", OutputVariableInfo::Element, OutputVariableInfo::Real, n_strain_components, _strain);
    femData->outController.setOutput(var1.name, var1);
    OutputVariableInfo var2("STRESS", OutputVariableInfo::Element, OutputVariableInfo::Real, n_strain_components, _stress);
    femData->outController.setOutput(var2.name, var2);
};

int FunctionStressStrain::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Solving ELEMENT_VELOCITY...");

    const MeshLib::IMesh *msh = _dis->getMesh();
    FemLib::FemNodalFunctionVector* u = (FemLib::FemNodalFunctionVector*)getInput(Displacement);
    FemLib::FEMIntegrationPointFunctionVector* strain = _strain;
    FemLib::FEMIntegrationPointFunctionVector* stress = _stress;
    FemLib::LagrangianFeObjectContainer* feObjects = u->getFeObjectContainer();

    const size_t dim = msh->getDimension();
    const size_t n_strain_components = getNumberOfStrainComponents();

    Ogs6FemData* femData = Ogs6FemData::getInstance();

    //calculate strain, stress
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++)
    {
        // element setup
        MeshLib::IElement* e = msh->getElemenet(i_e);
        const size_t nnodes = e->getNumberOfNodes();
        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        strain->setNumberOfIntegationPoints(i_e, n_gp);
        stress->setNumberOfIntegationPoints(i_e, n_gp);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());

        MaterialLib::Solid *solidphase = femData->list_solid[e->getGroupID()];
        // set D
        NumLib::LocalMatrix matD = NumLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
        //matD *= .0;
        NumLib::LocalMatrix nv(1,1);
        NumLib::LocalMatrix E(1,1);
        solidphase->poisson_ratio->eval(e_pos, nv);
        solidphase->Youngs_modulus->eval(e_pos, E);
        double Lambda, G, K;
        MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
        MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);

        // local u
        NumLib::LocalVector local_u(dim*nnodes);
        for (size_t j=0; j<nnodes; j++) {
            const size_t node_id = e->getNodeID(j);
            for (size_t k=0; k<dim; k++)
                local_u[j*dim+k] = u->getValue(node_id)(k);
        }

        // for each integration points
        NumLib::LocalMatrix matB(n_strain_components, nnodes*dim);
        NumLib::LocalMatrix matN(dim, nnodes*dim);
        double r[3] = {};
        double x[3] = {};
        for (size_t ip=0; ip<n_gp; ip++) {
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            NumLib::LocalMatrix &N = *fe->getBasisFunction();
            const NumLib::LocalMatrix &dN = *fe->getGradBasisFunction();
            fe->getRealCoordinates(x);

            // set N,B
            setNu_Matrix(dim, nnodes, N, matN);
            setB_Matrix(dim, nnodes, dN, matB);

            // strain
            NumLib::LocalVector gp_strain(n_strain_components);
            gp_strain.noalias() = matB * local_u;
            strain->setIntegrationPointValue(i_e, ip, gp_strain);

            // stress
            NumLib::LocalVector gp_stress(n_strain_components);
            gp_stress = matD * gp_strain;
            stress->setIntegrationPointValue(i_e, ip, gp_stress);

//                std::cout << "strain=\n" << gp_strain << std::endl;
//                std::cout << "D=\n" << matD << std::endl;
//                std::cout << "stress=\n" << gp_stress << std::endl;
        }
    }

    setOutput(Strain, strain);
    setOutput(Stress, stress);

    return 0;
}

