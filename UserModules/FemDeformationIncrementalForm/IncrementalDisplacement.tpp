/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IncrementalDisplacement.tpp
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#include "logog.hpp"

#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "MaterialLib/PorousMedia.h"
#include "MaterialLib/Solid.h"
#include "PhysicsLib/FemLinearElasticTools.h"
#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

template<class T1, class T2>
typename FunctionIncrementalDisplacement<T1, T2>::MyVariable* FunctionIncrementalDisplacement<
        T1, T2>::getDisplacementComponent(MyVariable *u_x, MyVariable* u_y,
        MyVariable* u_z, const std::string &var_name)
{
    if (var_name.find("_X") != std::string::npos) {
        return u_x;
    } else if (var_name.find("_Y") != std::string::npos) {
        return u_y;
    } else {
        return u_z;
    }
}

template<class T1, class T2>
bool FunctionIncrementalDisplacement<T1, T2>::initialize(
        const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //--------------------------------------------------------------------------
    // set up mesh and FE objects
    //--------------------------------------------------------------------------
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    const FemLib::PolynomialOrder::type msh_order = FemLib::PolynomialOrder::Linear;
    if (msh_order != FemLib::PolynomialOrder::Linear) {
        INFO("->generating higher order mesh...");
        MeshLib::MeshGenerator::generateHigherOrderUnstrucuredMesh(*(MeshLib::UnstructuredMesh*)msh, 2);
        INFO("* mesh id %d: order=%d, nodes=%d, elements=%d", msh_id, msh->getMaxiumOrder(), msh->getNumberOfNodes(msh->getMaxiumOrder()), msh->getNumberOfElements());
    }
    const size_t dim = msh->getDimension();
    const MeshLib::CoordinateSystem msh_coord =
            msh->getGeometricProperty()->getCoordinateSystem();
    MyDiscreteSystem* dis =
            DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);


    //--------------------------------------------------------------------------
    // set up problem
    //--------------------------------------------------------------------------
    _problem = new MyProblemType(dis);
    _problem->setTimeSteppingFunction(*tim);

    //--------------------------------------------------------------------------
    // set up variables
    //--------------------------------------------------------------------------
    // definitions
    MyVariable* u_x = _problem->addVariable("u_x", msh_order);
    u_x->setFeObjectContainer(_feObjects);
    MyVariable* u_y = _problem->addVariable("u_y", msh_order);
    u_y->setFeObjectContainer(_feObjects);
    MyVariable* u_z = NULL;
    if (dim==3) {
        u_z = _problem->addVariable("u_z", msh_order);
        u_z->setFeObjectContainer(_feObjects);
    }
    FemVariableBuilder var_builder;
    var_builder.doit(this->getOutputParameterName(Displacement) + "_X1", option,
            msh, femData->geo, femData->geo_unique_name, _feObjects, u_x);
    var_builder.doit(this->getOutputParameterName(Displacement) + "_Y1", option,
            msh, femData->geo, femData->geo_unique_name, _feObjects, u_y);
    if (dim==3) {
        var_builder.doit(this->getOutputParameterName(Displacement)+"_Z1", option, msh, femData->geo, femData->geo_unique_name, _feObjects, u_z);
    }

    //--------------------------------------------------------------------------
    // set up equations
    //--------------------------------------------------------------------------
    MyEquationType* eqs = _problem->createEquation();
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(
            _feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(
            _feObjects, msh_coord);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);
    eqs->initialize(linear_assembler, r_assembler, j_eqs);

    //--------------------------------------------------------------------------
    // set up this solution algorithm
    //--------------------------------------------------------------------------
    _solution = new MySolutionType(dis, _problem);
    typename MySolutionType::LinearSolverType* linear_solver =
            _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);
    DiscreteLib::DofEquationIdTable* dofMapping =
            _solution->getDofEquationIdTable();
    dofMapping->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);
    dofMapping->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_POINT);

    //--------------------------------------------------------------------------
    // set up stress, strain
    //--------------------------------------------------------------------------
    // create strain, stress vectors
    const size_t n_strain_components = getNumberOfStrainComponents(dim);
    MathLib::LocalVector v0 = MathLib::LocalVector::Zero(n_strain_components);
    _strain = new MyIntegrationPointFunctionVector();
    _strain->initialize(dis, v0);
    _previous_stress = new MyIntegrationPointFunctionVector();
    _previous_stress->initialize(dis, v0);
    _current_stress = new MyIntegrationPointFunctionVector(*_previous_stress);


    //--------------------------------------------------------------------------
    // set up data for output
    //--------------------------------------------------------------------------
    // create u variable which is vector
    MathLib::LocalVector tmp_u0(3);
    tmp_u0 *= .0;
    _displacement = new MyNodalFunctionVector();
    _displacement->initialize(*dis, u_x->getCurrentOrder(), tmp_u0);
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        for (size_t j=0; j<dim; j++) {
            _displacement->getValue(i)(j) = _solution->getCurrentSolution(j)->getValue(i);
        }
    }
    for (size_t i=0; i<dim; i++) {
        _vec_u_components.push_back(new NodalPointScalarWrapper(_displacement, i));
    }
    for (size_t i=0; i<n_strain_components; i++) {
        _vec_strain_components.push_back(new IntegrationPointScalarWrapper(_strain, i));
        _vec_stress_components.push_back(new IntegrationPointScalarWrapper(_current_stress, i));
    }

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id,
            OutputVariableInfo::Node, OutputVariableInfo::Real, 2,
            _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i = 0; i < _vec_u_components.size(); i++) {
        OutputVariableInfo var1(
                this->getOutputParameterName(Displacement)
                        + getDisplacementComponentPostfix(i), msh_id,
                OutputVariableInfo::Node, OutputVariableInfo::Real, 1,
                _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (size_t i=0; i<n_strain_components; i++) {
        OutputVariableInfo var1(this->getOutputParameterName(Strain) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_strain_components[i]);
        femData->outController.setOutput(var1.name, var1);
        OutputVariableInfo var2(this->getOutputParameterName(Stress) + getStressStrainComponentPostfix(i), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 1, _vec_stress_components[i]);
        femData->outController.setOutput(var2.name, var2);
    }

    //--------------------------------------------------------------------------
    // set initial output
    //--------------------------------------------------------------------------
    this->setOutput(Displacement, _displacement);
    this->setOutput(Strain, _strain);
    this->setOutput(Stress, _current_stress);

    return true;
}

template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::initializeTimeStep(
        const NumLib::TimeStep &/*time*/)
{
    INFO("->pre-processing in solving time step...");
    this->_problem->getEquation()->getResidualAssembler()->setPreviousStress(_previous_stress);
}

template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::updateOutputParameter(
        const NumLib::TimeStep &/*time*/)
{
    INFO("->post-processing in solving time step...");
    MeshLib::IMesh* msh = _problem->getDiscreteSystem()->getMesh();
    const size_t dim = msh->getDimension();
    //--------------------------------------------------------------------------
    // compute strain, stress
    //--------------------------------------------------------------------------
    for (size_t i=0; i<_displacement->getNumberOfNodes(); i++) {
        for (size_t j=0; j<dim; j++) {
            _displacement->getValue(i)(j) = _solution->getCurrentSolution(j)->getValue(i);
        }
    }
    calculateStrainStress(_displacement, _strain, _current_stress);

    //--------------------------------------------------------------------------
    // set system output
    //--------------------------------------------------------------------------
    setOutput(Displacement, _displacement);
    setOutput(Strain, _strain);
    setOutput(Stress, _current_stress);
}


template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::output(
        const NumLib::TimeStep &/*time*/)
{
    const size_t msh_id = _problem->getDiscreteSystem()->getMesh()->getID();
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Displacement), msh_id,
            OutputVariableInfo::Node, OutputVariableInfo::Real, 2,
            _displacement);
    femData->outController.setOutput(var.name, var);
    for (size_t i = 0; i < _vec_u_components.size(); i++) {
        OutputVariableInfo var1(
                this->getOutputParameterName(Displacement)
                        + getDisplacementComponentPostfix(i), msh_id,
                OutputVariableInfo::Node, OutputVariableInfo::Real, 1,
                _vec_u_components[i]);
        femData->outController.setOutput(var1.name, var1);
    }
}


template<class T1, class T2>
void FunctionIncrementalDisplacement<T1, T2>::calculateStrainStress(MyNodalFunctionVector* u, MyIntegrationPointFunctionVector* strain, MyIntegrationPointFunctionVector* stress)
{
    INFO("->computing strain and stress...");

    const MeshLib::IMesh* msh = _problem->getDiscreteSystem()->getMesh();
    FemLib::LagrangeFeObjectContainer* feObjects = _feObjects;

    const size_t dim = msh->getDimension();
    const size_t n_strain_components = getNumberOfStrainComponents(dim);

    Ogs6FemData* femData = Ogs6FemData::getInstance();

    //calculate strain, stress
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++)
    {
        // element setup
        MeshLib::IElement* e = msh->getElement(i_e);
        const size_t nnodes = e->getNumberOfNodes();
        FemLib::IFiniteElement *fe = feObjects->getFeObject(*e);
        FemLib::IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        strain->setNumberOfIntegationPoints(i_e, n_gp);
        stress->setNumberOfIntegationPoints(i_e, n_gp);
        const NumLib::TXPosition e_pos(NumLib::TXPosition::Element, e->getID());

        MaterialLib::Solid *solidphase = femData->list_solid[e->getGroupID()];
        // set D
        MathLib::LocalMatrix matD = MathLib::LocalMatrix::Zero(n_strain_components, n_strain_components);
        //matD *= .0;
        MathLib::LocalMatrix nv(1,1);
        MathLib::LocalMatrix E(1,1);
        solidphase->poisson_ratio->eval(e_pos, nv);
        solidphase->Youngs_modulus->eval(e_pos, E);
        double Lambda, G, K;
        MaterialLib::calculateLameConstant(nv(0,0), E(0,0), Lambda, G, K);
        MaterialLib::setElasticConsitutiveTensor(dim, Lambda, G, matD);

        // local u
        MathLib::LocalVector local_u(dim*nnodes);
        for (size_t j=0; j<nnodes; j++) {
            const size_t node_id = e->getNodeID(j);
            for (size_t k=0; k<dim; k++)
                local_u[j*dim+k] = u->getValue(node_id)(k);
        }

        // for each integration points
        MathLib::LocalMatrix matB = MathLib::LocalMatrix::Zero(n_strain_components, nnodes*dim);
        MathLib::LocalMatrix matN = MathLib::LocalMatrix::Zero(dim, nnodes*dim);
        double r[3] = {};
        double x[3] = {};
        for (size_t ip=0; ip<n_gp; ip++) {
            integral->getSamplingPoint(ip, r);
            fe->computeBasisFunctions(r);
            MathLib::LocalMatrix &N = *fe->getBasisFunction();
            const MathLib::LocalMatrix &dN = *fe->getGradBasisFunction();
            fe->getRealCoordinates(x);

            // set N,B
            setNu_Matrix_byPoint(dim, nnodes, N, matN);
            setB_Matrix_byPoint(dim, nnodes, dN, matB);

            // strain
            MathLib::LocalVector gp_strain(n_strain_components);
            gp_strain.noalias() = matB * local_u;
            strain->setIntegrationPointValue(i_e, ip, gp_strain);

            // stress
            MathLib::LocalVector gp_stress(n_strain_components);
            gp_stress = matD * gp_strain;
            stress->setIntegrationPointValue(i_e, ip, gp_stress);

//                std::cout << "strain=\n" << gp_strain << std::endl;
//                std::cout << "D=\n" << matD << std::endl;
//                std::cout << "stress=\n" << gp_stress << std::endl;
        }
    }
}


