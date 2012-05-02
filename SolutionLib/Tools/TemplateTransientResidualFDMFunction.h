
#pragma once

#include <vector>

#include "MathLib/Function/IFunction.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/TransientAssembler/StencilWiseTransientResidualAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FdmLib/FdmFunction.h"
#include "FdmLib/BoundaryConditions.h"


namespace SolutionLib
{

typedef DiscreteLib::DiscreteVector<double> MyFemVector;

/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_USER_FEM_PROBLEM
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
template <
	class T_USER_FEM_PROBLEM,
    class T_LOCAL_RESIDUAL_ASSEMBLER
    >
class TemplateTransientResidualFDMFunction
	: public MathLib::TemplateFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_USER_FEM_PROBLEM UserFemProblem;
    typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;

    /// constructor
    /// @param problem		Fem problem
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientResidualFDMFunction(UserFemProblem* problem, DiscreteLib::DofEquationIdTable* dofManager, UserLocalResidualAssembler* asssembler)
        : _problem(problem), _local_assembler(asssembler), _dofManager(dofManager),
          _t_n1(0), _u_n0(0), _st(0)
    {
    };

    ///
	virtual ~TemplateTransientResidualFDMFunction() {};

	///
    MathLib::TemplateFunction<MyFemVector,MyFemVector>* clone() const
	{
    	return new TemplateTransientResidualFDMFunction
    				<
						UserFemProblem,
						UserLocalResidualAssembler
    				>(_problem, _dofManager, _local_assembler);
	}

    /// reset property
    void reset(const NumLib::TimeStep* t, MyFemVector* u_n0, MyFemVector* st)
    {
        this->_t_n1 = const_cast<NumLib::TimeStep*>(t);
        this->_u_n0 = u_n0;
        this->_st = st;
    };

    /// evaluate residual
    /// @param u_n1	current results
    /// @param r residual
    void eval(const MyFemVector &u_n1, MyFemVector &r)
    {
    	// input, output
        const NumLib::TimeStep &t_n1 = *this->_t_n1;
        MyFemVector *u_n = this->_u_n0;

        // prepare data
        //UserFemProblem* pro = _problem;

        //TODO temporally
        std::vector<MyFemVector*> vec_un;
        vec_un.push_back(const_cast<MyFemVector*>(u_n));
        std::vector<MyFemVector*> vec_un1;
        vec_un1.push_back(const_cast<MyFemVector*>(&u_n1));

		// assembly
        NumLib::StencilWiseTransientResidualAssembler assembler(&t_n1, &vec_un, &vec_un1, _local_assembler);
        r = .0;
        assembler.assembly(*_problem->getMesh(), *_dofManager, r);

        // add source/sink term (i.e. RHS in linear equation)
        if (_st!=0)
        	r += *_st;
    }

private:
    UserFemProblem* _problem;
    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    NumLib::TimeStep* _t_n1;
    MyFemVector* _u_n0;
    MyFemVector* _st;
};


} //end
