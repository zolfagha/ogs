
#pragma once

#include <vector>

#include "MathLib/Function/IFunction.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/TransientAssembler/ElementWiseTransientResidualAssembler.h"
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
    class T_LOCAL_RESIDUAL_ASSEMBLER
    >
class TemplateTransientResidualFEMFunction
	: public MathLib::TemplateFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_USER_FEM_PROBLEM UserFemProblem;
    typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;

    /// constructor
    /// @param problem		Fem problem
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientResidualFEMFunction(UserFemProblem &problem, DiscreteLib::DofEquationIdTable &dofManager, UserLocalResidualAssembler &asssembler)
        : _problem(&problem), _local_assembler(&asssembler), _dofManager(&dofManager)
    {
    };

    ///
	virtual ~TemplateTransientResidualFEMFunction() {};

	///
    MathLib::TemplateFunction<MyFemVector,MyFemVector>* clone() const
	{
    	return new TemplateTransientResidualFEMFunction<
						UserFemProblem,
						UserLocalResidualAssembler
    				>(*_problem, *_dofManager, *_local_assembler);
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
        MyFemVector *u_n = 0;

        // prepare data
        UserFemProblem* pro = _problem;

        //TODO temporally
        std::vector<MyFemVector*> vec_un;
        vec_un.push_back(const_cast<MyFemVector*>(u_n));
        std::vector<MyFemVector*> vec_un1;
        vec_un1.push_back(const_cast<MyFemVector*>(&u_n1));

		// assembly
        NumLib::ElementWiseTransientResidualAssembler assembler(t_n1, vec_un, vec_un1, *_local_assembler);
        assembler.assembly(*_problem->getMesh(), *_dofManager, r);
    }


private:
    UserFemProblem* _problem;
    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    NumLib::TimeStep* _t_n1;
};


} //end
