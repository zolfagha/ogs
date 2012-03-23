
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/LinearEquation/MeshBasedDiscreteLinearEquation.h"
#include "DiscreteLib/Utils/SparsityBuilder.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "NumLib/TransientAssembler/DiscreteLinearEquationAssembler.h"

#include "ISolution.h"
#include "TransientFemFunctions.h"
#include "Nonlinear.h"

namespace SolutionLib
{


template<
	template <class> class T_TIME_ODE_ASSEMBLER,
    class T_LINEAR_SOLVER,
    class T_FEM_PROBLEM_AND_ASSEMBLER,
    class T_NONLINEAR
    >
class SingleStepFEM;


/**
 * \brief Solution algorithm for linear transient problems using FEM with Euler time stepping method
 *
 * @tparam T_TIME_ODE_ASSEMBLER Time ODE equation assembler
 * @tparam T_LINEAR_SOLVER     	Linear equation solver class
 * @tparam T_USER_FEM_PROBLEM  	FEM problem class
 * @tparam T_USER_ASSEMBLY     	Local assembly class
 */
template <
	template <class> class T_TIME_ODE_ASSEMBLER,
    class T_LINEAR_SOLVER,
    template <class> class T_USER_FEM_PROBLEM,
    class T_USER_ASSEMBLY,
    class T_NONLINEAR
    >
class SingleStepFEM<T_TIME_ODE_ASSEMBLER, T_LINEAR_SOLVER, T_USER_FEM_PROBLEM<T_USER_ASSEMBLY>, T_NONLINEAR>
	: public AbstractTimeSteppingAlgorithm
{
public:
    typedef T_TIME_ODE_ASSEMBLER<T_USER_ASSEMBLY> UserTimeOdeAssembler;
    typedef T_USER_FEM_PROBLEM<T_USER_ASSEMBLY> UserFemProblem;
    typedef TemplateTransientLinearFEMFunction<T_TIME_ODE_ASSEMBLER,T_LINEAR_SOLVER,T_USER_FEM_PROBLEM,T_USER_ASSEMBLY> UserLinearFemFunction;

    /// constructor
    ///
    /// - initialize solution vectors
    /// - set up DoFs and equation index
    /// - prepare linear equation and solver
    /// - prepare linear functions
    SingleStepFEM(DiscreteLib::DiscreteSystem &dis, UserFemProblem &problem)
        : AbstractTimeSteppingAlgorithm(*problem.getTimeSteppingFunction()), 
          _problem(&problem), _element_ode_assembler(problem.getElementAssemlby()),
          _discrete_system(&dis), _linear_eqs(0)
    {
        const size_t n_var = problem.getNumberOfVariables();
        // create dof map
        for (size_t i=0; i<n_var; i++) {
            _dofManager.addVariableDoFs(dis.getMesh()->getID(), 0, problem.getField(i)->getNumberOfNodes());
        }
        _dofManager.construct();
        // initialize solution vectors
//        _u_n.resize(n_var, 0);
        _u_n1.resize(n_var, 0);
        _u_n1[0] = (FemLib::FemNodalFunctionScalar*) problem.getIC(0)->clone();
//        for (size_t i=0; i<n_var; i++) {
//            _u_n1[i] = (FemLib::FemNodalFunctionScalar*) problem.getIC(i)->clone();
//        }
        size_t n_dofs = _dofManager.getTotalNumberOfActiveDoFs();
        _vec_n0 = dis.createVector<double>(n_dofs);
        _vec_n1 = dis.createVector<double>(_dofManager.getTotalNumberOfActiveDoFs());
        FemLib::FemNodalFunctionScalar *f_ic = (FemLib::FemNodalFunctionScalar*) problem.getIC(0); //TODO one var
        *_vec_n1 = *f_ic->getNodalValuesAsStdVec();
        // create linear equation systems
        _linear_solver = new T_LINEAR_SOLVER();
        _linear_eqs = _discrete_system->createLinearEquation<DiscreteLib::TemplateMeshBasedDiscreteLinearEquation, T_LINEAR_SOLVER, DiscreteLib::SparsityBuilderFromNodeConnectivity>(*_linear_solver, _dofManager);
        // setup linear function
        _linear_fucntion = new UserLinearFemFunction(problem, *_linear_eqs);
        //MathLib::NewtonFunctionDXVector<UserLinearFemFunction, T_LINEAR_SOLVER, MathLib::CRSMatrix<double, signed> > f_dx(*_linear_fucntion, _linear_solver);
        _nonlinear = new T_NONLINEAR(*_linear_fucntion);
    };

    ///
    virtual ~SingleStepFEM()
    {
        Base::releaseObject(_linear_solver);
        Base::releaseObject(_linear_fucntion);
        Base::releaseObject(_nonlinear);
    }

    /// get a linear equation solver object
    T_LINEAR_SOLVER* getLinearEquationSolver() { return _linear_solver; }

    /// solve 
    int solveTimeStep(const NumLib::TimeStep &t_n1)
    {
        double dt = t_n1.getTime() - AbstractTimeSteppingAlgorithm::getTimeStepFunction()->getPrevious();
        NumLib::TimeStep this_t_n1;
        this_t_n1.assign(t_n1);
        this_t_n1.setTimeStepSize(dt);
        *_vec_n0 = *_vec_n1;

        _linear_fucntion->reset(t_n1);
        _nonlinear->solve(*_linear_fucntion, *_vec_n0, *_vec_n1);

        _u_n1[0]->setNodalValues(*_vec_n1);

        return 0;
    }

    /// get the current solution
    FemLib::FemNodalFunctionScalar* getCurrentSolution(int var_id)
    {
        return _u_n1[var_id];
    }

    /// get the time ode assembler
    UserTimeOdeAssembler* getTimeODEAssembler() { return &_element_ode_assembler; }

private:
    UserFemProblem* _problem;
    UserTimeOdeAssembler _element_ode_assembler;
    DiscreteLib::DofEquationIdTable _dofManager;
    T_LINEAR_SOLVER* _linear_solver;
    DiscreteLib::DiscreteSystem *_discrete_system;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
//    std::vector<FemLib::FemNodalFunctionScalar*> _u_n;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n1;
    UserLinearFemFunction* _linear_fucntion;
    T_NONLINEAR* _nonlinear;
    MyFemVector *_vec_n0;
    MyFemVector *_vec_n1;

};


}
