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
#include "logog.hpp"

#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "Ogs6FemData.h"

#include "../FemDeformationTotalForm/FemLinearElasticTools.h"

void exactSol_Mode1(double xx, double yy, double k1, double kappa, double mu, double lambda, double uu, double vv)
{
    const static double pi = 4*atan(1.0);

    double rr = sqrt(xx*xx+yy*yy); // Crack tip is at (0, 0).
    double drrdx = xx/rr;
    double drrdy = yy/rr;

    double dthdx = -yy/(xx*xx+yy*yy);
    double dthdy = 0;
    if (xx == 0.)
        dthdy = 0;
    else
        dthdy = 1/(xx+yy*yy/xx);

    double th = 0;
    if (xx == 0) {
        if (yy == 0) {
            INFO("Theta is undetermined for (x,y) = (0,0)!");
            th = 0;
        } else if (yy > 0) {
            th = 0.5*pi;
        } else {
            th = 1.5*pi;
        }
    } else if (yy == 0) {
        if (xx > 0)
            th = 0;
        else
            th = pi;
    } else if (xx > 0) {
        th = atan(yy/xx);
    } else if (xx < 0) {
        th = pi+atan(yy/xx);
    } else {
        INFO("Internal error!");
    }

    if (th > pi) // This is important due to the multiplication with 1/2 below!!!
        th = th-2*pi;

    // Get exact solution.
    uu = k1/(2*mu)*sqrt(rr/(2*pi))*cos(0.5*th)*(kappa-1+2*sin(0.5*th)*sin(0.5*th));
    vv = k1/(2*mu)*sqrt(rr/(2*pi))*sin(0.5*th)*(kappa+1-2*cos(0.5*th)*cos(0.5*th));

//    duudrr = 1/8*k1*2^(1/2)*cos(1/2*th)*(kappa+1-2*cos(1/2*th)^2)/mu/pi^(1/2)/rr^(1/2);
//    duudth = -1/8*k1*2^(1/2)/pi^(1/2)*rr^(1/2)*sin(1/2*th)*(kappa+1-6*cos(1/2*th)^2)/mu;
//    duudx = duudrr * drrdx + duudth * dthdx;
//    duudy = duudrr * drrdy + duudth * dthdy;
//
//    dvvdrr = 1/8*k1/mu*2^(1/2)/pi^(1/2)/rr^(1/2)*sin(1/2*th)*(kappa+1-2*cos(1/2*th)^2);
//    dvvdth = 1/8*k1*2^(1/2)/pi^(1/2)*rr^(1/2)*cos(1/2*th)*(kappa+5-6*cos(1/2*th)^2)/mu;
//    dvvdx = dvvdrr * drrdx + dvvdth * dthdx;
//    dvvdy = dvvdrr * drrdy + dvvdth * dthdy;
//
//    Sigma11 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(1-sin(0.5*th)*sin(1.5*th));
//    Sigma22 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(1+sin(0.5*th)*sin(1.5*th));
//    Sigma12 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(sin(0.5*th)*cos(1.5*th));
//
//    % Sigma = [Sigma11 Sigma12; Sigma12 Sigma22];
//    % Eps = -0.25*lambda/(mu*(lambda+mu))*trace(Sigma)*[1 0; 0 1] + 1/(2*mu)*Sigma;
//    % Eps11 = Eps(1,1);
//    % Eps12 = Eps(1,2);
//    % Eps22 = Eps(2,2);
//
//    Eps11 = duudx;
//    Eps12 = 0.5 * (duudy + dvvdx);
//    Eps22 = dvvdy;
}

void searchMinMaxNodes(const MeshLib::IMesh &msh, std::vector<size_t> &found_nodes)
{
    double x_min=1e+99, x_max = -1e+99;
    double y_min=1e+99, y_max = -1e+99;
    double pt[3];
    //search x min/max
    for (size_t i=0; i<msh.getNumberOfNodes(); i++) {
      const GeoLib::Point *pt = msh.getNodeCoordinatesRef(i);
      if ((*pt)[0]<x_min) x_min = (*pt)[0];
      if ((*pt)[0]>x_max) x_max = (*pt)[0];
      if ((*pt)[1]<y_min) y_min = (*pt)[1];
      if ((*pt)[1]>y_max) y_max = (*pt)[1];
    }

    //search nodes on min/max
    for (size_t i=0; i<msh.getNumberOfNodes(); i++) {
        const GeoLib::Point *pt = msh.getNodeCoordinatesRef(i);
      if (abs((*pt)[0]-x_min)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      } else if (abs((*pt)[0]-x_max)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      } else if (abs((*pt)[1]-y_min)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      } else if (abs((*pt)[1]-y_max)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      }
    }
}

void getEnrichedNodesElems(const MeshLib::IMesh &msh, const std::vector<double> &ff, std::vector<size_t> &ElemsEnriched, std::vector<size_t> &NodesEnriched)
{
    // Get cut elements and correspoding nodes.

    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        MeshLib::IElement* e = msh.getElemenet(i);
        bool hasNegativeNode = false;
        bool hasNonNegativeNode = false;
        size_t cnt = 0;
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (ff[e->getNodeID(j)] < .0) hasNegativeNode = true;
            else hasNonNegativeNode = true;
            if (msh.getNodeCoordinatesRef(e->getNodeID(j))->getData()[0]<.0) {
                cnt++;
            }
        }
        if (hasNegativeNode && hasNonNegativeNode && cnt==e->getNumberOfNodes())
            ElemsEnriched.push_back(i);
    }

    std::set<size_t> set_nodes;
    for (size_t i=0; i<ElemsEnriched.size(); i++) {
        MeshLib::IElement* e = msh.getElemenet(ElemsEnriched[i]);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            set_nodes.insert(e->getNodeID(j));
        }
    }
    NodesEnriched.assign(set_nodes.begin(), set_nodes.end());
}

template <class T1, class T2>
bool FunctionDisplacement<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    _msh = femData->list_mesh[msh_id];
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
    double k1 = 1.0;
    double EE = 10000.;
    double nu = 0.3;
    double fx = 0, fy = 0;
    double mu = EE/(2.*(1.+nu));
    //lambda = EE*nu/((1+nu)*(1-2*nu)); % plane strain
    //kappa = 3-4*nu; % plane strain
    double lambda = EE*nu/(1.-nu*nu); // plane stress
    double kappa  = (3.-nu)/(1.+nu); // plane stress

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
    std::vector<double> uDirValues;
    std::vector<double> vDirValues;
    for (size_t i=0; i<Bound.size(); i++) {
        uDirValues[i] = uuExact[Bound[i]];
        vDirValues[i] = vvExact[Bound[i]];
    }

    // Get enriched elements and nodes.
    std::vector<size_t> ElemsEnriched, NodesEnriched;
    getEnrichedNodesElems(*_msh, ff, ElemsEnriched, NodesEnriched);

    // initialize LinearEQS
    MathLib::DenseLinearEquation leqs;
    leqs.create(4*NodeNum);

    // domain integration
    for (size_t i=0; i<ElemNum; i++) {

    }

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

