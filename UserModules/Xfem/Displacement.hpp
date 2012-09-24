/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Displacement.hpp
 *
 * Created on 2012-09-20 by Norihiro Watanabe
 */

#include <cmath>
#include <limits>
#include <vector>

#include "logog.hpp"

#include "MeshLib/Tools/MeshGenerator.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "Ogs6FemData.h"

#include "../FemDeformationTotalForm/FemLinearElasticTools.h"
#include "matlab_functions.h"

namespace xfem
{

template <class T1, class T2>
bool FunctionDisplacement<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
//    size_t msh_id = option.getOption<size_t>("MeshID");
//    size_t time_id = option.getOption<size_t>("TimeGroupID");
//    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    _msh = MeshLib::MeshGenerator::generateRegularQuadMesh(2., 19, -1., -1., .0);
    _dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(_msh);
    _feObjects = new FemLib::LagrangianFeObjectContainer(*_msh);

#if 0
    // set up problem
    _problem = new MyProblemType(_dis);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    MyVariable* u_x = _problem->addVariable("u_x");
    MyVariable* u_y = _problem->addVariable("u_y");
    // IC
    NumLib::TXFunctionBuilder f_builder;
    MyNodalFunctionScalar* u0 = new MyNodalFunctionScalar();
    u0->initialize(*_dis, FemLib::PolynomialOrder::Linear, 0);
    u_x->setIC(u0);
    u_y->setIC(u0);
    // BC
    const BaseLib::Options* opBCList = option.getSubGroup("BCList");
    for (const BaseLib::Options* opBC = opBCList->getFirstSubGroup("BC"); opBC!=0; opBC = opBCList->getNextSubGroup())
    {
        std::string var_name = opBC->getOption("Variable");
        std::string geo_type = opBC->getOption("GeometryType");
        std::string geo_name = opBC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string dis_name = opBC->getOption("DistributionType");
        double dis_v = opBC->getOption<double>("DistributionValue");
        NumLib::ITXFunction* f_bc =  f_builder.create(dis_name, dis_v);
        getDisplacementComponent(u_x, u_y, 0, var_name)->addDirichletBC(new SolutionLib::FemDirichletBC(_msh, geo_obj, f_bc));
    }

    // ST
    const BaseLib::Options* opSTList = option.getSubGroup("STList");
    for (const BaseLib::Options* opST = opSTList->getFirstSubGroup("ST"); opST!=0; opST = opSTList->getNextSubGroup())
    {
        std::string var_name = opST->getOption("Variable");
        std::string geo_type = opST->getOption("GeometryType");
        std::string geo_name = opST->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string st_type = opST->getOption("STType");
        std::string dis_name = opST->getOption("DistributionType");
        double dis_v = opST->getOption<double>("DistributionValue");
        if (st_type.compare("NEUMANN")==0) {
            dis_v *= -1; //TODO
        }
        NumLib::ITXFunction* f_st =  f_builder.create(dis_name, dis_v);
        if (f_st!=NULL) {
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(_msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(_msh, geo_obj, f_st);
            }
            getDisplacementComponent(u_x, u_y, 0, var_name)->addNeumannBC(femSt);
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }
#endif

//    // set initial output
//    OutputVariableInfo var(this->getOutputParameterName(Displacement), OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
//    femData->outController.setOutput(var.name, var);
//    for (size_t i=0; i<_vec_u_components.size(); i++) {
//        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
//
//    // initial output parameter
//    this->setOutput(Displacement, _displacement);


    return true;
}


template <class T1, class T2>
int FunctionDisplacement<T1,T2>::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    const size_t NodeNum = _msh->getNumberOfNodes();
    const size_t ElemNum = _msh->getNumberOfElements();

    // define test case parameters
    const double k1 = 1.0;
    const double EE = 10000.;
    const double nu = 0.3;
    const double fx = 0, fy = 0;
    const double mu = EE/(2.*(1.+nu));
    //lambda = EE*nu/((1+nu)*(1-2*nu)); % plane strain
    //kappa = 3-4*nu; % plane strain
    const double lambda = EE*nu/(1.-nu*nu); // plane stress
    const double kappa  = (3.-nu)/(1.+nu); // plane stress

    // get level-set function
    std::vector<double> ff(NodeNum);
    for (size_t i=0; i<NodeNum; i++)
        ff[i] = _msh->getNodeCoordinatesRef(i)->getData()[1];

    // Get exact solution at the nodes.
    std::vector<double> uuExact(NodeNum, .0), vvExact(NodeNum, .0);
    for (size_t i=0; i<NodeNum; i++) {
        const GeoLib::Point *pt = _msh->getNodeCoordinatesRef(i);
        exactSol_Mode1((*pt)[0], (*pt)[1], k1, kappa, mu, lambda, uuExact[i], vvExact[i]);
    }

    // define Dirichlet BC
    std::vector<size_t> Bound;
    searchMinMaxNodes(*_msh, Bound);
    std::vector<size_t> uDirNodes(Bound), vDirNodes(Bound);
    for (size_t i=0; i<vDirNodes.size(); i++)
        vDirNodes[i] += NodeNum;
    std::vector<double> uDirValues(Bound.size()), vDirValues(Bound.size());
    for (size_t i=0; i<Bound.size(); i++) {
        uDirValues[i] = uuExact[Bound[i]];
        vDirValues[i] = vvExact[Bound[i]];
    }

    // Get enriched elements and nodes.
    std::vector<size_t> ElemsEnriched, NodesEnriched;
    getEnrichedNodesElems(*_msh, ff, ElemsEnriched, NodesEnriched);
    std::set<size_t> SetNodesEnriched;
    SetNodesEnriched.insert(NodesEnriched.begin(), NodesEnriched.end());

    // initialize LinearEQS
    MathLib::DenseLinearEquation leqs;
    leqs.create(4*NodeNum);

    // domain integration
    for (size_t i=0; i<ElemNum; i++) {

        MeshLib::IElement* e = _msh->getElement(i);
        const size_t n_ele_nodes = e->getNumberOfNodes();
        NumLib::LocalVector Nodes(n_ele_nodes);
        FemLib::IFiniteElement* fe = _feObjects->getFeObject(*e);
        NumLib::LocalVector xxElem(n_ele_nodes), yyElem(n_ele_nodes);
        NumLib::LocalVector ffEle(n_ele_nodes);
        for (size_t j=0; j<n_ele_nodes; j++) {
            Nodes(j) = e->getNodeID(j);
            xxElem(j) = _msh->getNodeCoordinatesRef(e->getNodeID(j))->getData()[0];
            yyElem(j) = _msh->getNodeCoordinatesRef(e->getNodeID(j))->getData()[1];
            ffEle(j) = ff[e->getNodeID(j)];
        }

        // activate nodes are enriched
        NumLib::LocalVector NodesAct(n_ele_nodes);
        for (size_t i=0; i<n_ele_nodes; i++) {
            NodesAct(i) = (SetNodesEnriched.count(Nodes(i))>0 ? 1 : 0);
        }

        // set integration points in the reference element
        FemLib::IFemNumericalIntegration *q = fe->getIntegrationMethod();
        const size_t nQnQ = q->getNumberOfSamplingPoints();
        std::vector<GeoLib::Point> vec_int_ref_xx(nQnQ);
        std::vector<double> vec_int_ref_w(nQnQ);
        for (size_t j=0; j<nQnQ; j++) {
            q->getSamplingPoint(j, vec_int_ref_xx[j].getData());
            vec_int_ref_w[j] = q->getWeight(j);
        }
        NumLib::LocalVector xxIntRef, yyIntRef, wwIntRef;
        IntPoints2DLevelSet(ffEle, vec_int_ref_xx, vec_int_ref_w, xxIntRef, yyIntRef, wwIntRef);
        const size_t Curr_nQ = xxIntRef.rows();

        // get shape functions
        NumLib::LocalMatrix N, dNdx, dNdy;
        NumLib::LocalMatrix M, dMdx, dMdy;
        NumLib::LocalVector xxInt, yyInt, wwInt, ffInt;
        ShapeFctsXFEMSign(
                xxElem, yyElem, ffEle, NodesAct, xxIntRef, yyIntRef, wwIntRef,
                Curr_nQ,
                N, dNdx, dNdy, M, dMdx, dMdy, xxInt, yyInt, wwInt, ffInt);

        // integrate
        BuildMatRhs_Hooke(
                N, dNdx, dNdy, M, dMdx, dMdy,
                xxInt, yyInt, wwInt, ffInt, Nodes,
                lambda, lambda, mu, mu, fx, fy, Curr_nQ, NodeNum,
                leqs);
    }

    // Insert Dirichlet BCs.
    leqs.setKnownX(uDirNodes, uDirValues);
    leqs.setKnownX(vDirNodes, vDirValues);

//    // Reduce system of equations.
//    Pos = [[1:1:2*NodeNum]'; NodesEnriched+2*NodeNum; NodesEnriched+3*NodeNum];
//    MAT = MAT(Pos, Pos);
//    RHS = RHS(Pos);

    //disp(sprintf('Condition number of final system         : %15.5e', condest(MAT)))

    // Solve system of equations for solution.
    leqs.solve();
    double *x = leqs.getX();

    leqs.printout(std::cout);
//    uuTotal = Sol(1:NodeNum);
//    vvTotal = Sol(NodeNum+1:2*NodeNum);

    return 0;
}

template <class T1, class T2>
void FunctionDisplacement<T1,T2>::accept(const NumLib::TimeStep &/*time*/)
{
//    //update data for output
//    const size_t n_strain_components = getNumberOfStrainComponents();
//    Ogs6FemData* femData = Ogs6FemData::getInstance();
//    for (size_t i=0; i<n_strain_components; i++) {
//        OutputVariableInfo var1(this->getOutputParameterName(NodStrain) + getStressStrainComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_strain_components[i]);
//        femData->outController.setOutput(var1.name, var1);
//        OutputVariableInfo var2(this->getOutputParameterName(NodStress) + getStressStrainComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_nodal_stress_components[i]);
//        femData->outController.setOutput(var2.name, var2);
//    }
}

template <class T1, class T2>
void FunctionDisplacement<T1,T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
//    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
//        _displacement->getValue(i)(0) = _solution->getCurrentSolution(0)->getValue(i);
//        _displacement->getValue(i)(1) = _solution->getCurrentSolution(1)->getValue(i);
//    }
    setOutput(Displacement, _displacement);

    //calculateStressStrain();
}

template <class T1, class T2>
void FunctionDisplacement<T1,T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement), OutputVariableInfo::Node, OutputVariableInfo::Real, 2, _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i=0; i<_vec_u_components.size(); i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Displacement) + getDisplacementComponentPostfix(i), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
};

}
