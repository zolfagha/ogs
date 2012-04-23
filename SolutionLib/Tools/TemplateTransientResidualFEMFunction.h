
#pragma once

#include <vector>

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
 * \brief Template class for FEM residual functions
 *
 * \tparam T_USER_FEM_PROBLEM
 * \tparam T_LOCAL_ASSEMBLER
 * \tparam T_LINEAR_SOLVER
 */
template <
	class T_USER_FEM_PROBLEM,
    class T_LOCAL_LINEAR_ASSEMBLER,
    class T_LINEAR_SOLVER
    >
class TemplateTransientResidualFEMFunction
	: public MathLib::TemplateFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_USER_FEM_PROBLEM UserFemProblem;
    typedef T_LOCAL_LINEAR_ASSEMBLER UserLocalLinearAssembler;

    /// constructor
    /// @param problem		Fem problem
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientResidualFEMFunction(UserFemProblem &problem, UserLocalLinearAssembler &asssembler, DiscreteLib::IDiscreteLinearEquation &linear_eqs)
        : _problem(&problem), _local_assembler(&asssembler),  _linear_eqs(&linear_eqs)
    {
    };

    ///
	virtual ~TemplateTransientResidualFEMFunction() {};

	///
    MathLib::TemplateFunction<MyFemVector,MyFemVector>* clone() const
	{
    	return new TemplateTransientLinearFEMFunction<
						UserFemProblem,
						UserLocalLinearAssembler,
						T_LINEAR_SOLVER
    				>(*_problem, *_local_assembler, *_linear_eqs);
	}

    /// reset property
    void reset(const NumLib::TimeStep &t)
    {
    	this->_t_n1 = const_cast<NumLib::TimeStep*>(&t);
    };

    /// evaluate residual
    /// @param u_n1	current results
    /// @param r residual
    void eval(const MyFemVector &u_n1, MyFemVector &r)
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
        vec_un1.push_back(const_cast<MyFemVector*>(&u_n1));

		// assembly
        _linear_eqs->initialize();
        NumLib::ElementBasedTransientAssembler assembler(t_n1, vec_un, vec_un1, *_local_assembler);
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
    UserLocalLinearAssembler *_local_assembler;
    DiscreteLib::IDiscreteLinearEquation* _linear_eqs;
    NumLib::TimeStep* _t_n1;
};


} //end
