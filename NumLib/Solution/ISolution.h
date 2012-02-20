
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/DiscreteSystem.h"
#include "NumLib/Solution/IProblem.h"


namespace NumLib
{

/**
 * \brief Interface to all solution algorithm classes
 */
class ISolutionAlgorithm
{
public:
};

/**
 * \brief Abstract class for time-stepping method
 */
class AbstractTimeSteppingAlgorithm : public ISolutionAlgorithm, public ITransientSystem
{
private:
    ITimeStepFunction* _tim;

public:
    AbstractTimeSteppingAlgorithm(ITimeStepFunction &tim) : _tim(&tim) {};

    ITimeStepFunction* getTimeStepFunction() const {return _tim;};

    double suggestNext(const TimeStep &time_current)
    {
        return _tim->next(time_current.getTime());
    }

    bool isAwake(const TimeStep &time)
    {
        return (time.getTime()==_tim->next(time.getTime()));
    }

};

/**
 * \brief Solution algorithm for linear transient problems using FEM with Euler time stepping method
 */
template<class T_USER_FEM_PROBLEM, class T_USER_ASSEMBLY>
class TimeEulerSpaceFemLinearAlgorithm : public AbstractTimeSteppingAlgorithm
{
public:
    ///
    TimeEulerSpaceFemLinearAlgorithm(T_USER_FEM_PROBLEM &problem) : AbstractTimeSteppingAlgorithm(*problem.getTimeSteppingFunction()), _problem(&problem), _element_ode_assembler(problem.getElementAssemlby())
    {
        size_t n_var = problem.getNumberOfVariables();
        _u_n1.resize(n_var, 0);
        _u_n.resize(n_var, 0);
        for (size_t i=0; i<n_var; i++) {
            _u_n[i] = (FemLib::FemNodalFunctionScalar*) problem.getIC(i);
        }
    };

    /// solve 
    int solveTimeStep(const TimeStep &t_n1)
    {
        return solve(t_n1, _u_n, _u_n1);
    }

    FemLib::FemNodalFunctionScalar* getCurrentSolution(int var_id)
    {
        return _u_n1[var_id];
    }

private:
    /// 
	int solve(const TimeStep &t_n1, const std::vector<FemLib::FemNodalFunctionScalar*>& u_n, std::vector<FemLib::FemNodalFunctionScalar*>& u_n1) 
    {
        FemIVBVProblem<T_USER_ASSEMBLY> *pro = _problem;
        size_t n_var = pro->getNumberOfVariables();
        for (size_t i=0; i<n_var; i++) {
            u_n1[i] = new FemLib::FemNodalFunctionScalar(*u_n[i]);
        }
		
		// collect data
		double delta_t = t_n1.getTimeStep();
        double theta = 1.0;
		

        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) 
            pro->getFemDirichletBC(i)->setup();
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) 
            pro->getFemNeumannBC(i)->setup();

		// assembly
        const MeshLib::IMesh *msh = u_n1[0]->getMesh();
        IDiscreteLinearEquation* discreteEqs = _discrete_system->getLinearEquation();
        discreteEqs->construct(ElementBasedAssembler(TimeEulerElementAssembler<T_USER_ASSEMBLY>(t_n1, theta, _element_ode_assembler)));

        //apply BC1,2
        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) 
            pro->getFemDirichletBC(i)->apply(*discreteEqs->getLinearEquation());
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) 
            pro->getFemNeumannBC(i)->apply(discreteEqs->getRHS());
		
		// solve
		discreteEqs->solve();
        double *x = discreteEqs->getX();

        for (size_t i=0; i<n_var; i++) {
            u_n1[i]->setNodalValues(x); //TODO use DOF
        }

		// check 

        return 0;
	}

private:
    T_USER_FEM_PROBLEM* _problem;
    T_USER_ASSEMBLY _element_ode_assembler;

    DiscreteSystem *_discrete_system;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n;
    std::vector<FemLib::FemNodalFunctionScalar*> _u_n1;
};

}
