
#pragma once

#include <vector>

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TransientAssembler/ElementWiseTransientResidualAssembler.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemNeumannBC.h"

#include "FemVariable.h"

namespace SolutionLib
{

typedef DiscreteLib::DiscreteVector<double> MyFemVector;

/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
template <
    class T_LOCAL_RESIDUAL_ASSEMBLER
    >
class TemplateTransientResidualFEMFunction
	: public NumLib::TemplateFunction<MyFemVector, MyFemVector>
{
public:
    typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;

    /// constructor
    /// @param problem		Fem problem
    /// @param linear_eqs	Discrete linear equation
    TemplateTransientResidualFEMFunction(MeshLib::IMesh* msh, std::vector<FemVariable*>* list_var, DiscreteLib::DofEquationIdTable* dofManager, UserLocalResidualAssembler* asssembler)
        : _local_assembler(asssembler), _dofManager(dofManager),
          _t_n1(0), _u_n0(0), _st(0), _msh(msh)
    {
    };

    ///
	virtual ~TemplateTransientResidualFEMFunction() {};

	///
	NumLib::TemplateFunction<MyFemVector,MyFemVector>* clone() const
	{
    	return new TemplateTransientResidualFEMFunction
    				<
						UserLocalResidualAssembler
    				>(_list_var, _dofManager, _local_assembler);
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

        //TODO temporally
        std::vector<MyFemVector*> vec_un;
        vec_un.push_back(const_cast<MyFemVector*>(u_n));
        std::vector<MyFemVector*> vec_un1;
        vec_un1.push_back(const_cast<MyFemVector*>(&u_n1));

		// assembly
        NumLib::ElementWiseTransientResidualAssembler assembler(&t_n1, &vec_un, &vec_un1, _local_assembler);
        r = .0;
        assembler.assembly(*_msh, *_dofManager, r);

        // add source/sink term (i.e. RHS in linear equation)
        if (_st!=0)
        	r += *_st;
    }

private:
    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    NumLib::TimeStep* _t_n1;
    MyFemVector* _u_n0;
    MyFemVector* _st;
    std::vector<FemVariable*>* _list_var;
    MeshLib::IMesh* _msh;
};


} //end
