/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepFEM.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "DiscreteLib/Serial/DiscreteLinearEquation.h"
#include "DiscreteLib/Utils/SparsityBuilder.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TransientAssembler/ElementWiseTransientLinearEQSAssembler.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/DataType.h"
#include "SolutionLib/FemProblem/FemEquation.h"
#include "SolutionLib/FemProblem/FemDirichletBC.h"
#include "SolutionLib/FemProblem/FemNeumannBC.h"

#include "AbstractTimeSteppingAlgorithm.h"

namespace SolutionLib
{

/**
 * \brief Solution algorithm for linear transient problems using FEM with single time stepping method
 *
 * @tparam T_USER_FEM_PROBLEM      FEM problem class
 * @tparam T_LINEAR_SOLVER         Linear equation solver class
 */
template <
    class T_USER_FEM_PROBLEM,
    class T_LINEAR_SOLVER
    >
class SingleStepFEM
    : public AbstractTimeSteppingAlgorithm
{
public:
    typedef T_USER_FEM_PROBLEM UserFemProblem;
    typedef typename UserFemProblem::EquationType::LinearEQSType UserLinearFunction;
    typedef typename UserFemProblem::EquationType::ResidualEQSType UserResidualFunction;
    typedef typename UserFemProblem::EquationType::DxEQSType UserDxFunction;
    typedef NumLib::TemplateDiscreteNonlinearSolver<UserLinearFunction, UserResidualFunction, UserDxFunction> NonlinearSolverType;
    typedef T_LINEAR_SOLVER LinearSolverType;
    //typedef DiscreteLib::DiscreteVector<double> ImplVector;

    /// constructor
    /// - initialize solution vectors
    /// - set up DoFs and equation index
    /// - prepare linear equation and solver
    /// - prepare linear functions
    SingleStepFEM(DiscreteLib::DiscreteSystem* dis, UserFemProblem* problem)
        : AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()), 
          _problem(problem), _discrete_system(dis)
    {
        const size_t n_var = problem->getNumberOfVariables();
        MeshLib::IMesh* msh = dis->getMesh();

        // create dof map
        for (size_t i=0; i<n_var; i++) {
            size_t n_dof_per_var = msh->getNumberOfNodes();
            _dofManager.addVariableDoFs(msh->getID(), 0, n_dof_per_var);
        }
        _dofManager.construct(DiscreteLib::DofNumberingType::BY_POINT);
        const size_t n_total_dofs = _dofManager.getTotalNumberOfActiveDoFs();

        // initialize vectors for each variable
        _vec_u_n1.resize(n_var, 0);
        for (size_t i=0; i<n_var; i++) {
            _vec_u_n1[i] = problem->getVariable(i)->getIC()->clone();
        }

        // initialize vectors for solution, ST
        _x_n0 = dis->createVector<double>(n_total_dofs);
        _x_n1 = dis->createVector<double>(n_total_dofs);
        _x_n1_0 = dis->createVector<double>(n_total_dofs);
        _x_st = dis->createVector<double>(n_total_dofs);

        // copy values of each variable to one solution vector
        for (size_t i=0; i<n_var; i++) {
            SolutionVector* vec_var = problem->getVariable(i)->getIC()->getNodalValues();
            DiscreteLib::setGlobalVector(_dofManager, i, msh->getID(), *vec_var, *_x_n1);
        }

        // create linear equation systems
        _linear_solver = new LinearSolverType();
        _linear_eqs = _discrete_system->createLinearEquation
                <   LinearSolverType,
                    DiscreteLib::SparsityBuilderFromNodeConnectivity
                >(*_linear_solver, _dofManager);

        // setup functions
        std::vector<FemVariable*> list_var(n_var);
        for (size_t i=0; i<n_var; i++) list_var[i] = problem->getVariable(i);
        _f_linear = new UserLinearFunction(list_var, problem->getEquation()->getLinearAssembler(), _linear_eqs);
        _f_r = new UserResidualFunction(dis->getMesh(), list_var, &_dofManager, problem->getEquation()->getResidualAssembler());
        _f_dx = new UserDxFunction(dis->getMesh(), list_var, problem->getEquation()->getJacobianAssembler(), _linear_eqs);

        // setup nonlinear solver
        _f_nonlinear = new NonlinearSolverType(dis, _f_linear, _f_r, _f_dx);
    };

    ///
    virtual ~SingleStepFEM()
    {
        _discrete_system->deleteLinearEquation(_linear_eqs);
        _discrete_system->deleteVector(_x_n0);
        _discrete_system->deleteVector(_x_n1);
        _discrete_system->deleteVector(_x_n1_0);
        _discrete_system->deleteVector(_x_st);
        BaseLib::releaseObject(_linear_solver);
        BaseLib::releaseObject(_f_linear);
        BaseLib::releaseObject(_f_r);
        BaseLib::releaseObject(_f_dx);
        BaseLib::releaseObject(_f_nonlinear);
    }

    /// solve 
    int solveTimeStep(const NumLib::TimeStep &t_n1)
    {
        // time step
        double dt = t_n1.getTime() - AbstractTimeSteppingAlgorithm::getTimeStepFunction()->getPrevious();
        NumLib::TimeStep this_t_n1;
        this_t_n1.assign(t_n1);
        this_t_n1.setTimeStepSize(dt);

        const size_t n_var = _problem->getNumberOfVariables();
        const size_t msh_id = _discrete_system->getMesh()->getID();
        // bc1
        std::vector<size_t> list_bc1_eqs_id;
        std::vector<double> list_bc1_val;
        for (size_t i_var=0; i_var<n_var; i_var++) {
            FemVariable* var = _problem->getVariable(i_var);
            for (size_t i=0; i<var->getNumberOfDirichletBC(); i++) {
                SolutionLib::FemDirichletBC *bc1 = var->getDirichletBC(i);
                bc1->setup();
                std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
                std::vector<double> &list_bc_values = bc1->getListOfBCValues();
                DiscreteLib::convertToEqsValues(_dofManager, i_var, msh_id, list_bc_nodes, list_bc_values, list_bc1_eqs_id, list_bc1_val);
            }
        }

        // st
        std::vector<size_t> list_st_eqs_id;
        std::vector<double> list_st_val;
        for (size_t i_var=0; i_var<n_var; i_var++) {
            FemVariable* var = _problem->getVariable(i_var);
            for (size_t i=0; i<var->getNumberOfNeumannBC(); i++) {
                SolutionLib::IFemNeumannBC *bc2 = var->getNeumannBC(i);
                bc2->setup();
                std::vector<size_t> &list_bc_nodes = bc2->getListOfBCNodes();
                std::vector<double> &list_bc_values = bc2->getListOfBCValues();

                const size_t msh_id = _problem->getMesh()->getID();
                DiscreteLib::convertToEqsValues(_dofManager, i_var, msh_id, list_bc_nodes, list_bc_values, list_st_eqs_id, list_st_val);
            }
        }
        (*_x_st) = .0;
        for (size_t i=0; i<list_st_eqs_id.size(); i++) {
            (*_x_st)[list_st_eqs_id[i]] = list_st_val[i];
        }

        // setup functions
        _f_linear->reset(&t_n1, _x_n0);
        _f_r->reset(&t_n1, _x_n0, _x_st);
        _f_dx->reset(&t_n1, _x_n0);

        // initial guess
        *_x_n1_0 = *_x_n0;
        for (size_t i=0; i<list_bc1_eqs_id.size(); i++) {
            (*_x_n1_0)[list_bc1_eqs_id[i]] = list_bc1_val[i];
        }
        
        // solve
        _f_nonlinear->solve(*_x_n1_0, *_x_n1);

        // distribute solution vector to local vector for each variable
        for (size_t i=0; i<n_var; i++) {
            //SolutionVector* vec_var = _problem->getVariable(i)->getIC()->getNodalValues();
            DiscreteLib::setLocalVector(_dofManager, i, msh_id, *_x_n1, *_vec_u_n1[i]->getNodalValues());
        }

        return 0;
    }

    /// get the current solution
    FemLib::FemNodalFunctionScalar* getCurrentSolution(size_t var_id) { return _vec_u_n1[var_id]; }

    ///
    UserFemProblem* getProblem() {return _problem;};

    /// get a linear solver
    LinearSolverType* getLinearEquationSolver() { return _linear_solver; };
    /// get a nonlinear solver
    NonlinearSolverType* getNonlinearSolver() { return _f_nonlinear;};

    ///
    virtual void accept(const NumLib::TimeStep &t)
    {
        AbstractTimeSteppingAlgorithm::accept(t);
        *_x_n0 = *_x_n1; //copy current value to previous value
    };

private:
    DISALLOW_COPY_AND_ASSIGN(SingleStepFEM);

    //void convertToEqsValues(const DiscreteLib::DofEquationIdTable &eqs_map, size_t var_id, size_t msh_id, const std::vector<size_t> &list_node_id, const std::vector<double> &list_node_values, std::vector<size_t> &list_eqs_id, std::vector<double> &list_eqs_val)
    //{
    //    for (size_t j=0; j<list_node_id.size(); j++) {
    //        size_t pt_id = list_node_id[j];
    //        if (eqs_map.isActiveDoF(var_id, msh_id, pt_id)) {
    //            size_t eqs_id = eqs_map.mapEqsID(var_id, msh_id, pt_id);
    //            double bc_val = list_node_values[j];
    //            list_eqs_id.push_back(eqs_id);
    //            list_eqs_val.push_back(bc_val);
    //        }
    //    }
    //}


private:
    DiscreteLib::DofEquationIdTable _dofManager;
    UserFemProblem* _problem;
    LinearSolverType* _linear_solver;
    DiscreteLib::DiscreteSystem *_discrete_system;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    std::vector<FemLib::FemNodalFunctionScalar*> _vec_u_n1;
    UserLinearFunction* _f_linear;
    UserResidualFunction* _f_r;
    UserDxFunction* _f_dx;
    NonlinearSolverType* _f_nonlinear;
    SolutionVector *_x_n0;
    SolutionVector *_x_n1_0;
    SolutionVector *_x_n1;
    SolutionVector *_x_st;
};


}
