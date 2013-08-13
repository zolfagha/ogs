/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemResidualEQS.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Utils/Tools.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TransientAssembler/TransientElementWiseVectorUpdater.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "FemLib/Function/FemFunction.h"

#include "SolutionLib/Fem/FemVariable.h"
#include "SolutionLib/DataType.h"

#include "UserModules/FemReactGIAReduct/ReductConc.h"


/**
 * \brief Template class for FEM residual functions
 *
 * \tparam T_LOCAL_RESIDUAL_ASSEMBLER
 */
template <class T_DIS_SYS, class T_USER_FUNCTION_DATA>
class TemplateTransientResidualFEMFunction_GIA_Reduct
    : public NumLib::TemplateFunction<SolutionLib::SolutionVector, SolutionLib::SolutionVector>
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef SolutionLib::FemVariable MyVariable;
//    typedef T_LOCAL_RESIDUAL_ASSEMBLER UserLocalResidualAssembler;
//    typedef typename NumLib::TransientElementWiseVectorUpdater<UserLocalResidualAssembler> MyUpdater;
//    typedef typename T_DIS_SYS::template MyVectorAssembler<double, MyUpdater>::type MyGlobalAssembler;


    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
    TemplateTransientResidualFEMFunction_GIA_Reduct(MyDiscreteSystem* sys,
    		                                        const std::vector<MyVariable*> &list_var,
    		                                        DiscreteLib::DofEquationIdTable* dofManager,
    		                                        T_USER_FUNCTION_DATA* function_data)
        : _dis_sys(sys), _dofManager(dofManager),
          _t_n1(0), _u_n0(0), _st(0), _list_var(list_var), _function_data(function_data)
    {
    };

    ///
    virtual ~TemplateTransientResidualFEMFunction_GIA_Reduct() {};

    ///
    NumLib::TemplateFunction<SolutionLib::SolutionVector,SolutionLib::SolutionVector>* clone() const
    {
        return new TemplateTransientResidualFEMFunction_GIA_Reduct
                    <T_DIS_SYS, T_USER_FUNCTION_DATA>(_dis_sys, _list_var, _dofManager, _function_data);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const SolutionLib::SolutionVector* u_n0, SolutionLib::SolutionVector* st)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
        this->_st = st;
    };

    /// evaluate residual
    /// @param u_n1    current results
    /// @param r residual
    void eval(const SolutionLib::SolutionVector &u_n1, SolutionLib::SolutionVector &r);

	const NumLib::TimeStep* getTimeStepObj(void)
	{
		return _t_n1; 	
	}


private:
    MyDiscreteSystem *_dis_sys;
//    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    const NumLib::TimeStep* _t_n1;
    const SolutionLib::SolutionVector* _u_n0;
    SolutionLib::SolutionVector* _st;
    std::vector<MyVariable*> _list_var;
    T_USER_FUNCTION_DATA* _function_data;
};


template <class T_DIS_SYS, class T_USER_FUNCTION_DATA>
void TemplateTransientResidualFEMFunction_GIA_Reduct<T_DIS_SYS, T_USER_FUNCTION_DATA>::eval(const SolutionLib::SolutionVector &u_n1, SolutionLib::SolutionVector &r)
{
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
    const SolutionLib::SolutionVector *u_n = this->_u_n0;
    size_t msh_id = _dis_sys->getMesh()->getID();
    size_t node_idx, i;


    // assembly
    MeshLib::IMesh* msh = _dis_sys->getMesh();
    r = .0;
    //node based operations


    //FunctionReductConc *test = template<msh, ggf> FunctionReductConc();

	// calculate the global residual
    _function_data->GlobalResidualAssembler(t_n1, u_n1, r);
    //FunctionReductConc::GlobalResidualAssembler(&t_n1, &u_n1, &r);


//    // element based operation for assemblyCh
//    MyUpdater updater(&t_n1, msh, _dofManager, u_n, &u_n1, _local_assembler);
//    MyGlobalAssembler assembler(&updater);
//    r.construct(*_dofManager, assembler);



    // add source/sink term (i.e. RHS in linear equation)
    if (_st!=0) {
        r += *_st;
    }

    // set residuals to zero for Dirichlet BC
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
        	SolutionLib::FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            DiscreteLib::convertToEqsValues(*_dofManager, i, msh_id, bc1->getListOfBCNodes(), bc1->getListOfBCValues(), list_bc1_eqs_id, list_bc1_val);
        }
    }
    r.setSubvector(list_bc1_eqs_id, .0);
}
