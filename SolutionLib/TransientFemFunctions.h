
#pragma once

#include "MathLib/Function/IFunction.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"


namespace SolutionLib
{

typedef DiscreteLib::DiscreteVector<double> MyFemVector;

/**
 * \brief Template class for transient linear FEM functions
 */
template <
	template <class> class T_USER_FEM_PROBLEM,
	template <class> class T_TIME_ODE_ASSEMBLER,
    class T_LINEAR_SOLVER,
    class T_USER_ASSEMBLY
    >
class TemplateTransientLinearFEMFunction
	: public MathLib::IFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_TIME_ODE_ASSEMBLER<T_USER_ASSEMBLY> UserTimeOdeAssembler;
    typedef T_USER_FEM_PROBLEM<T_USER_ASSEMBLY> UserFemProblem;

    /// constructor
    /// @param problem		Fem problem
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientLinearFEMFunction(UserFemProblem &problem, DiscreteLib::IDiscreteLinearEquation &linear_eqs)
        : _problem(&problem), _element_ode_assembler(problem.getElementAssemlby()),
          _linear_eqs(&linear_eqs)
    {
    };

    ///
	virtual ~TemplateTransientLinearFEMFunction() {};

	///
    MathLib::IFunction<MyFemVector,MyFemVector>* clone() const
	{
    	return new TemplateTransientLinearFEMFunction<T_USER_FEM_PROBLEM, T_TIME_ODE_ASSEMBLER,T_LINEAR_SOLVER,T_USER_ASSEMBLY>(*_problem, *_linear_eqs);
	}

    /// reset property
    void reset(const NumLib::TimeStep &t)
    {
    	this->_t_n1 = const_cast<NumLib::TimeStep*>(&t);
    };

    /// solve linear equations discretized with FEM
    /// @param u_n	previous results
    /// @param u_n1 new results
    void eval(const MyFemVector &u_n, MyFemVector &u_n1)
    {
    	// input, output
        const NumLib::TimeStep &t_n1 = *this->_t_n1;

        // prepare data
        UserFemProblem* pro = _problem;

        // setup BC
        for (size_t i=0; i<pro->getNumberOfDirichletBC(); i++) {
            FemLib::FemDirichletBC<double> *bc1 = pro->getFemDirichletBC(i);
            bc1->setup();
            size_t varid = 0; //?
            _linear_eqs->setPrescribedDoF(varid, bc1->getListOfBCNodes(), bc1->getListOfBCValues());
        }
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++)
            pro->getFemNeumannBC(i)->setup();

        //TODO temporally
        std::vector<MyFemVector*> vec_un;
        vec_un.push_back(const_cast<MyFemVector*>(&u_n));
        std::vector<MyFemVector*> vec_un1;
        vec_un.push_back(const_cast<MyFemVector*>(&u_n1));

		// assembly
        _linear_eqs->initialize();
        NumLib::ElementBasedTransientAssembler assembler(t_n1, vec_un, vec_un1, _element_ode_assembler);
        _linear_eqs->construct(assembler);

        //apply BC1,2
        for (size_t i=0; i<pro->getNumberOfNeumannBC(); i++) {
            FemLib::FemNeumannBC<double, double> *bc2 = pro->getFemNeumannBC(i);
            size_t varid = 0; //?
            _linear_eqs->addRHS(varid, bc2->getListOfBCNodes(), bc2->getListOfBCValues(), -1.0);
        }

		// solve
		_linear_eqs->solve();
        _linear_eqs->getX(u_n1);
    }


private:
    UserFemProblem* _problem;
    UserTimeOdeAssembler _element_ode_assembler;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    NumLib::TimeStep* _t_n1;
};


} //end
