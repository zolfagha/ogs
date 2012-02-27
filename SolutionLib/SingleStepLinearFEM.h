
#pragma once

#include "ISolution.h"

#include <vector>

#include "Base/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MeshLib/Core/IMesh.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "DiscreteLib/DoF.h"
#include "DiscreteLib/DiscreteSystem.h"
#include "DiscreteLib/SparsityBuilder.h"

#include "NumLib/Discrete/DiscreteLinearEquationAssembler.h"

namespace SolutionLib
{

template<   template <class> class T_TIME_ODE_ASSEMBLER, 
            class T_LINEAR_SOLVER, 
            class T_FEM_PROBLEM_AND_ASSEMBLER >
class SingleStepLinearFEM;


/**
 * \brief Solution algorithm for linear transient problems using FEM with Euler time stepping method
 *
 * @tparam T_LINEAR_SOLVER     Linear equation solver class
 * @tparam T_USER_FEM_PROBLEM  FEM problem class
 * @tparam T_USER_ASSEMBLY     Local assembly class
 */
template<   template <class> class T_TIME_ODE_ASSEMBLER, 
            class T_LINEAR_SOLVER, 
            template <class> class T_USER_FEM_PROBLEM, 
            class T_USER_ASSEMBLY >
class SingleStepLinearFEM<T_TIME_ODE_ASSEMBLER, T_LINEAR_SOLVER, T_USER_FEM_PROBLEM<T_USER_ASSEMBLY>> : public AbstractTimeSteppingAlgorithm
{
public:
    typedef T_TIME_ODE_ASSEMBLER<T_USER_ASSEMBLY> UserTimeOdeAssembler;
    typedef T_USER_FEM_PROBLEM<T_USER_ASSEMBLY> UserFemProblem;

    ///
    SingleStepLinearFEM(DiscreteLib::DiscreteSystem &dis, UserFemProblem &problem) 
        : AbstractTimeSteppingAlgorithm(*problem.getTimeSteppingFunction()), 
        _discrete_system(&dis), _problem(&problem), _element_ode_assembler(problem.getElementAssemlby()), _linear_eqs(0)
    {
        const size_t n_var = problem.getNumberOfVariables();
        // initialize solution vectors
        _u_n.resize(n_var, 0);
        _u_n1.resize(n_var, 0);
        for (size_t i=0; i<n_var; i++) {
            _u_n1[i] = (FemLib::FemNodalFunctionScalar*) problem.getIC(i)->clone();
        }
        // create linear solver
        _linear_solver = new T_LINEAR_SOLVER();
        // create dof map
        for (size_t i=0; i<n_var; i++) {
            _dofManager.addDoF(problem.getField(i)->getNumberOfNodes());
        }
        _dofManager.construct();
        // create linear equation systems
        _linear_eqs = _discrete_system->createLinearEquation<T_LINEAR_SOLVER, DiscreteLib::SparsityBuilderFromNodeConnectivity>(*_linear_solver, _dofManager);
    };

    virtual ~SingleStepLinearFEM()
    {
        Base::releaseObject(_linear_solver);
    }

    T_LINEAR_SOLVER* getLinearEquationSolver()
    {
        return _linear_solver;
    }

    /// solve 
    int solveTimeStep(const NumLib::TimeStep &t_n1)
    {
        double dt = t_n1.getTime() - AbstractTimeSteppingAlgorithm::getTimeStepFunction()->getPrevious();
        NumLib::TimeStep this_t_n1;
        this_t_n1.assign(t_n1);
        this_t_n1.setTimeStepSize(dt);
        const size_t n_var = _problem->getNumberOfVariables();
        for (size_t i=0; i<n_var; i++) {
            if (_u_n[i]!=0) delete _u_n[i];
            _u_n[i] = _u_n1[i];
        }
        return solve(this_t_n1, _u_n, _u_n1);
    }

    /// get the current solution
    FemLib::FemNodalFunctionScalar* getCurrentSolution(int var_id)
    {
        return _u_n1[var_id];
    }

    /// get the time ode assembler
    UserTimeOdeAssembler* getTimeODEAssembler()
    {
        return &_element_ode_assembler;
    }

private:
    /// 
	int solve(const NumLib::TimeStep &t_n1, const std::vector<FemLib::FemNodalFunctionScalar*>& u_n, std::vector<FemLib::FemNodalFunctionScalar*>& u_n1) 
    {
        // prepare data
        UserFemProblem* pro = _problem;
        const size_t n_var = pro->getNumberOfVariables();
        for (size_t i=0; i<n_var; i++) {
            u_n1[i] = new FemLib::FemNodalFunctionScalar(*u_n[i]);
        }
        std::vector<DiscreteLib::DiscreteVector<double>*> vec_un(n_var);
        for (size_t i=0; i<n_var; i++) {
            vec_un[i] = u_n[i]->getNodalValuesAsStdVec();
        }

        // setup BC
        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) {
            FemLib::FemDirichletBC<double> *bc1 = pro->getFemDirichletBC(i);
            bc1->setup();
            size_t varid = 0; //?
            _linear_eqs->setPrescribedDoF(varid, bc1->getListOfBCNodes(), bc1->getListOfBCValues());
        }
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) 
            pro->getFemNeumannBC(i)->setup();

		// assembly
        _linear_eqs->construct(NumLib::ElementBasedTransientAssembler(t_n1, vec_un, _element_ode_assembler));

        //apply BC1,2
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) {
            FemLib::FemNeumannBC<double, double> *bc2 = pro->getFemNeumannBC(i);
            size_t varid = 0; //?
            _linear_eqs->addRHS(varid, bc2->getListOfBCNodes(), bc2->getListOfBCValues(), -1.0);
        }
		
		// solve
		_linear_eqs->solve();
        double *x = _linear_eqs->getLocalX();

        for (size_t i=0; i<n_var; i++) {
            u_n1[i]->setNodalValues(x); //TODO use DOF
        }

		// check 

        return 0;
	}

private:
    UserFemProblem* _problem;
    UserTimeOdeAssembler _element_ode_assembler;
    DiscreteLib::DofMapManager _dofManager;
    T_LINEAR_SOLVER* _linear_solver;
    DiscreteLib::DiscreteSystem *_discrete_system;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n1;
};


}
