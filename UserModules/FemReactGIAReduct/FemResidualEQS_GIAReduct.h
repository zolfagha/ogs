/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemResidualEQS.h
 *
 * Created on 2013-08-016 by Reza Zolfaghari & Norihiro Watanabe
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
#include "ReductionGIANodeInfo.h"


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
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
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
          _time_step(0), _u_n0(0), _st(0), _list_var(list_var), _function_data(function_data), _ReductionGIA(function_data->getReductionGIA()),
          _n_xi_Kin_bar(_ReductionGIA->get_n_xi_Kin_bar()), _n_xi_Mob(_ReductionGIA->get_n_xi_Mob()), _n_eta(_ReductionGIA->get_n_eta()), _n_eta_bar(_ReductionGIA->get_n_eta_bar()),
          _n_xi_Sorp_tilde(_ReductionGIA->get_n_xi_Sorp_tilde()), _n_xi_Min_tilde(_ReductionGIA->get_n_xi_Min_tilde()), _n_xi_Sorp(_ReductionGIA->get_n_xi_Sorp()), _n_xi_Min(_ReductionGIA->get_n_xi_Min()),
          _n_xi_Kin(_ReductionGIA->get_n_xi_Kin()),_n_xi_Min_bar(_ReductionGIA->get_n_xi_Min_bar()), _n_xi_local(_ReductionGIA->get_n_xi_local()), _n_xi_global(_ReductionGIA->get_n_xi_global()),
          _n_xi_Sorp_bar(_ReductionGIA->get_n_xi_Sorp_bar()), _J_tot_kin(_ReductionGIA->get_n_xi_Kin_total()), _n_xi_Sorp_bar_li(_ReductionGIA->get_n_xi_Sorp_bar_li()), _n_xi_Sorp_bar_ld(_ReductionGIA->get_n_xi_Sorp_bar_ld()),
          _xi_global(function_data->get_xi_global()), _xi_local(function_data->get_xi_local()), _eta(function_data->get_eta()), _eta_bar(function_data->get_eta_bar()),
          _global_vec_Rate(function_data->get_global_vec_Rate())
    {
    };

    ///
    virtual ~TemplateTransientResidualFEMFunction_GIA_Reduct() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    ///
    NumLib::TemplateFunction<SolutionLib::SolutionVector,SolutionLib::SolutionVector>* clone() const
    {
        return new TemplateTransientResidualFEMFunction_GIA_Reduct
                    <T_DIS_SYS, T_USER_FUNCTION_DATA>(_dis_sys, _list_var, _dofManager, _function_data);
    }

    /// reset property
    void reset(const NumLib::TimeStep& t, const SolutionLib::SolutionVector* u_n0, SolutionLib::SolutionVector* st)
    {
        _time_step = &t;
        _u_n0 = u_n0;
        _st = st;
    };

    /// evaluate residual
    /// @param u_n1    current results
    /// @param r residual
    void eval(const SolutionLib::SolutionVector &u_n1, SolutionLib::SolutionVector &r);

	const NumLib::TimeStep* getTimeStepObj(void)
	{
		return _time_step;
	}

private:
	void GlobalResidualAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, SolutionLib::SolutionVector & residual_global);
	void assembly(const NumLib::TimeStep & delta_t, const std::size_t _n_xi, const SolutionLib::SolutionVector & u_cur_xiglob,
	                                            SolutionLib::SolutionVector & residual_global, std::size_t switchOn, MathLib::LocalMatrix &_mat_Asorp, MathLib::LocalMatrix &_mat_Amin);

private:
    MyDiscreteSystem *_dis_sys;
//    UserLocalResidualAssembler *_local_assembler;
    DiscreteLib::DofEquationIdTable* _dofManager;
    const NumLib::TimeStep* _time_step;
    const SolutionLib::SolutionVector* _st;
    const SolutionLib::SolutionVector *_u_n0;
    std::vector<MyVariable*> _list_var;
    T_USER_FUNCTION_DATA* _function_data;

    //TODO pass via constructor
    ogsChem::chemReductionGIA* _ReductionGIA;
    // std::map<size_t, ReductionGIANodeInfo*>* _bc_info;
    std::vector<MyNodalFunctionScalar*> &_xi_global, &_xi_local, &_eta, &_eta_bar, &_global_vec_Rate;
    NumLib::ITXFunction* _vel;
    FemLib::IFemNumericalIntegration* _q;
    FemLib::IFiniteElement* _fe;
    std::size_t _n_xi_global, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Sorp, _n_xi_Min, _n_xi_Kin, _n_xi_local, _n_xi_Sorp_bar, _n_xi_Min_bar, _n_eta, _n_eta_bar, _n_xi_Mob, _n_xi_Kin_bar, _J_tot_kin
    			,_n_xi_Sorp_bar_li, _n_xi_Sorp_bar_ld;
};


template <class T_DIS_SYS, class T_USER_FUNCTION_DATA>
void TemplateTransientResidualFEMFunction_GIA_Reduct<T_DIS_SYS, T_USER_FUNCTION_DATA>::eval(const SolutionLib::SolutionVector &u_n1, SolutionLib::SolutionVector &r)
{
    // input, output
    const SolutionLib::SolutionVector *u_n = this->_u_n0;
    size_t msh_id = _dis_sys->getMesh()->getID();
    size_t node_idx, i;


    // assembly
    MeshLib::IMesh* msh = _dis_sys->getMesh();
    r = .0;
    //node based operations


    //FunctionReductConc *test = template<msh, ggf> FunctionReductConc();

	// calculate the global residual
    GlobalResidualAssembler(*this->_time_step, u_n1, r);


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

template <class T_DIS_SYS, class T_USER_FUNCTION_DATA>
void TemplateTransientResidualFEMFunction_GIA_Reduct<T_DIS_SYS, T_USER_FUNCTION_DATA>::GlobalResidualAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, SolutionLib::SolutionVector & residual_global)
{
    const size_t nnodes = _dis_sys->getMesh()->getNumberOfNodes();

    // current xi global
    MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde,
                         loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;
                         //global_cur_xi_Sorp_tilde, global_cur_xi_Min_tilde,
                         //global_cur_xi_Sorp, global_cur_xi_Min, global_cur_xi_Kin;
    // previous xi global
    MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;
                         //global_pre_xi_Sorp_tilde, global_pre_xi_Min_tilde, global_pre_xi_Sorp, global_pre_xi_Min, global_pre_xi_Kin;

    // current xi local
    MathLib::LocalVector loc_cur_xi_local, loc_xi_local, loc_cur_xi_Mob, loc_cur_xi_Sorp_bar, loc_cur_xi_Min_bar, loc_cur_xi_Kin_bar, loc_cur_xi_Sorp_bar_li, loc_cur_xi_Sorp_bar_ld;

    // residual vectors
    MathLib::LocalVector res43, res44, res45, res46, local_residual;//, residual_global_res43, residual_global_res44, residual_global_res45, residual_global_res46, global_vec_LHS_sorp, global_vec_RHS_sorp, global_vec_LHS_min, global_vec_RHS_min;

    // eta vectors
    MathLib::LocalVector loc_cur_eta, loc_cur_eta_bar;

    // rate vector
    MathLib::LocalVector vec_Rate, vec_Rate_45, vec_Rate_46; //, global_vec_Rate_45, global_vec_Rate_46;

    MathLib::LocalMatrix _mat_A1sorp, _mat_A2sorpli, _mat_A2sorpld, _mat_A1min, _mat_Ald;

    // initialize the local vector
    //current xi global
    loc_cur_xi_global               = MathLib::LocalVector::Zero( _n_xi_global );
    loc_cur_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde );
    loc_cur_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde );
    loc_cur_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp);
    loc_cur_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min);
    loc_cur_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin);
    //previous xi global
    loc_pre_xi_global               = MathLib::LocalVector::Zero( _n_xi_global );
    loc_pre_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde );
    loc_pre_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde );
    loc_pre_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp);
    loc_pre_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min);
    loc_pre_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin);
    //current xi local
    loc_cur_xi_local                = MathLib::LocalVector::Zero( _n_xi_local );
    loc_cur_xi_Sorp_bar             = MathLib::LocalVector::Zero( _n_xi_Sorp_bar );
    loc_cur_xi_Min_bar              = MathLib::LocalVector::Zero( _n_xi_Min_bar );
    loc_cur_xi_Sorp_bar_li          = MathLib::LocalVector::Zero(_n_xi_Sorp_bar_li );
    loc_cur_xi_Sorp_bar_ld          = MathLib::LocalVector::Zero(_n_xi_Sorp_bar_ld );
    //residual vectors
    res43                       = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde);
    res44                       = MathLib::LocalVector::Zero( _n_xi_Min_tilde);
    res45                       = MathLib::LocalVector::Zero( _n_xi_Sorp);
    res46                       = MathLib::LocalVector::Zero( _n_xi_Min);
//  residual_global_res43       = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde * nnodes);
//  residual_global_res44       = MathLib::LocalVector::Zero( _n_xi_Min_tilde * nnodes);
//  residual_global_res45       = MathLib::LocalVector::Zero( _n_xi_Min * nnodes);
//  residual_global_res46       = MathLib::LocalVector::Zero( _n_xi_Sorp * nnodes);
//  global_vec_LHS_sorp         = MathLib::LocalVector::Zero( _n_xi_Sorp * nnodes);
//  global_vec_RHS_sorp         = MathLib::LocalVector::Zero( _n_xi_Sorp * nnodes);
//  global_vec_LHS_min          = MathLib::LocalVector::Zero( _n_xi_Min * nnodes);
//  global_vec_RHS_min          = MathLib::LocalVector::Zero( _n_xi_Min * nnodes);
    // current eta mobie and immobile
    loc_cur_eta             = MathLib::LocalVector::Zero( _n_eta );
    loc_cur_eta_bar         = MathLib::LocalVector::Zero( _n_eta_bar );


    //initialize it before?
    _mat_A1sorp        = _ReductionGIA->get_matrix_A1sorp();
    _mat_A2sorpli      = _ReductionGIA->get_matrix_A2sorpli();
    _mat_A1min         = _ReductionGIA->get_matrix_A1min();
    _mat_Ald           = _ReductionGIA->get_matrix_Ald();
    _mat_A2sorpld      = _ReductionGIA->get_matrix_A2sorpld();

    MathLib::LocalMatrix _mat_Asorp = _mat_A1sorp - _mat_A2sorpli;
    MathLib::LocalMatrix _mat_Amin = _mat_A1min - _mat_Ald * _mat_A2sorpld;


    // loop over all the nodes
    for (size_t node_idx = 0; node_idx < nnodes; node_idx++ )
    {

            // on each node, get the right start value
            // get the right set of xis

            for (size_t i=0; i < _n_xi_global; i++)
                loc_cur_xi_global[i] = u_cur_xiglob[node_idx * _n_xi_global + i];
            for (size_t i=0; i < _n_xi_global; i++)
                loc_pre_xi_global[i] = _xi_global[i]->getValue(node_idx);
            for (size_t i=0; i < _n_xi_local; i++)
                loc_cur_xi_local[i] = _xi_local[i]->getValue(node_idx);
            for (size_t i=0; i < _n_eta; i++)
                loc_cur_eta[i] = _eta[i]->getValue(node_idx);
            // fill in eta_immob
            for (size_t i=0; i < _n_eta_bar; i++)
                loc_cur_eta_bar[i] = _eta_bar[i]->getValue(node_idx);


            // current xi global
            loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde);
            loc_cur_xi_Min_tilde    = loc_cur_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
            loc_cur_xi_Sorp         = loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
            loc_cur_xi_Min          = loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);
            loc_cur_xi_Kin          = loc_cur_xi_global.tail(_n_xi_Kin);


//          for (i=0; i < _n_xi_Sorp_tilde; i++)
//              _global_cur_xi_Sorp_tilde[i]->setValue(node_idx, loc_cur_xi_Sorp_tilde[i]);
//          for (i=0; i < _n_xi_Min_tilde; i++)
//              _global_cur_xi_Min_tilde[i]->setValue(node_idx, loc_cur_xi_Min_tilde[i]);
//          for (i=0; i < _n_xi_Sorp; i++)
//              _global_cur_xi_Sorp[i]->setValue(node_idx, loc_cur_xi_Sorp[i]);
//          for (i=0; i < _n_xi_Min; i++)
//              _global_cur_xi_Min[i]->setValue(node_idx, loc_cur_xi_Min[i]);
//          for (i=0; i < _n_xi_Kin; i++)
//              _global_cur_xi_Kin[i]->setValue(node_idx, loc_cur_xi_Kin[i]);


            // previous xi global
            loc_pre_xi_Sorp_tilde   = loc_pre_xi_global.head(_n_xi_Sorp_tilde);
            loc_pre_xi_Min_tilde    = loc_pre_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
            loc_pre_xi_Sorp         = loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
            loc_pre_xi_Min          = loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);
            loc_pre_xi_Kin          = loc_pre_xi_global.tail(_n_xi_Kin);

//          for (i=0; i < _n_xi_Sorp_tilde; i++)
//              _global_pre_xi_Sorp_tilde[i]->setValue(node_idx, loc_pre_xi_Sorp_tilde[i]);
//          for (i=0; i < _n_xi_Min_tilde; i++)
//              _global_pre_xi_Min_tilde[i]->setValue(node_idx, loc_pre_xi_Min_tilde[i]);
//          for (i=0; i < _n_xi_Sorp; i++)
//              _global_pre_xi_Sorp[i]->setValue(node_idx, loc_pre_xi_Sorp[i]);
//          for (i=0; i < _n_xi_Min; i++)
//              _global_pre_xi_Min[i]->setValue(node_idx, loc_pre_xi_Min[i]);
//          for (i=0; i < _n_xi_Kin; i++)
//              _global_pre_xi_Kin[i]->setValue(node_idx, loc_pre_xi_Kin[i]);


            // current xi local
            loc_cur_xi_Mob          = loc_cur_xi_local.head(_n_xi_Mob);
            loc_cur_xi_Sorp_bar     = loc_cur_xi_local.segment(_n_xi_Mob, _n_xi_Sorp_bar);
            loc_cur_xi_Min_bar      = loc_cur_xi_local.segment( _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar);
            loc_cur_xi_Kin_bar      = loc_cur_xi_local.tail(_n_xi_Kin_bar);

            loc_cur_xi_Sorp_bar_li  = loc_cur_xi_Sorp_bar.topRows(_n_xi_Sorp_bar_li);
            loc_cur_xi_Sorp_bar_ld  = loc_cur_xi_Sorp_bar.bottomRows(_n_xi_Sorp_bar_ld);





            // calculate node based AE (Eq. 3.43 and 3.44)
            res43 = loc_cur_xi_Sorp_tilde - loc_cur_xi_Sorp + loc_cur_xi_Sorp_bar_li;
            for(std::size_t i = 0; i < _n_xi_Sorp_tilde; i++)
            residual_global[_n_xi_global * node_idx + i] = res43[i];

            if(_n_xi_Sorp_bar_ld != 0)
            	res44 = loc_cur_xi_Min_tilde - loc_cur_xi_Min + loc_cur_xi_Min_bar + loc_cur_xi_Sorp_bar_ld;
            else
            	res44 = loc_cur_xi_Min_tilde - loc_cur_xi_Min + loc_cur_xi_Min_bar;

            for(std::size_t i = 0; i < _n_xi_Min_tilde; i++)
            residual_global[_n_xi_global * node_idx + i] = res44[i];

//          // collect the xi_local_new  ??
//          for (i=0; i < _n_xi_Sorp_tilde; i++)
//              residual_global_res43[i]->setValue(node_idx, res43[i]);
//          // residual_global_res43.push_back( res43 );
//
//          for (i=0; i < _n_xi_Min_tilde; i++)
//              residual_global_res44[i]->setValue(node_idx, res44[i]);
//          // residual_global_res44.push_back( res44 );


            // Note: since rate is already used in local problem, is it the same? probably not.
            // calculate the nodal kinetic reaction rates
            _ReductionGIA->Calc_Kin_Rate(loc_cur_xi_Mob,
                                            loc_pre_xi_Sorp,
                                            loc_cur_xi_Sorp_tilde,
                                            loc_cur_xi_Sorp_bar,
                                            loc_cur_xi_Min,
                                            loc_cur_xi_Min_tilde,
                                            loc_cur_xi_Min_bar,
                                            loc_cur_xi_Kin,
                                            loc_cur_xi_Kin_bar,
                                            loc_cur_eta,
                                            loc_cur_eta_bar,
                                            vec_Rate);



            //_J_tot_kin  = vec_Rate.rows();
            for (size_t i=0; i < _J_tot_kin; i++)
                _global_vec_Rate[i]->setValue(node_idx, vec_Rate[i]);

//          for (i=0; i < vec_Rate_46.cols(); i++)
//              _global_vec_Rate_46[i]->setValue(node_idx, vec_Rate_46[i]);


    }

    // solve for Eq. 45
    std::size_t switchOn = 1;  // 1 for xi sorp tilde
    this->assembly(delta_t, _n_xi_Sorp, u_cur_xiglob, residual_global, switchOn, _mat_Asorp, _mat_Amin);

    switchOn = 0;  // 1 for xi sorp tilde
    this->assembly(delta_t, _n_xi_Min, u_cur_xiglob, residual_global, switchOn, _mat_Asorp, _mat_Amin);

//    for(std::vector<int>::iterator iter; iter _residual_global_res45.begin(); iter != _residual_global_res45.end(); ++iter)
//    _residual_global_res45[iter] = _global_vec_LHS_sorp[iter] - _global_vec_RHS_sorp[iter] - _global_vec_Rate_45[iter];
//
//    // solve for Eq. 46
//    this->assembly(delta_t, _n_xi_Min, global_cur_xi_Min_tilde, global_pre_xi_Min_tilde, global_cur_xi_Min, global_pre_xi_Min, _global_vec_LHS_min, _global_vec_RHS_min);
//
//    residual_global_res46 = _global_vec_LHS_min - _global_vec_RHS_min - _global_vec_Rate_46;
//
//
//
//    // construct global residual vector
//  // loop over all the nodes
//    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
//       node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
//       node_idx++ )
//  {
//
//  for (i=0; i < _n_xi_Sorp_tilde; i++)
//      res43[i] = this->_residual_global_res43[i]->getValue(node_idx);
//  for (i=0; i < _n_xi_Min_tilde; i++)
//      res44[i] = this->_residual_global_res44[i]->getValue(node_idx);
//
//  for (i=0; i < _n_xi_Sorp; i++)
//      res45[i] = this->_residual_global_res45[i]->getValue(node_idx);
//  for (i=0; i < _n_xi_Min; i++)
//      res46[i] = this->_residual_global_res46[i]->getValue(node_idx);
//
//  local_residual.head(_n_xi_Sorp_tilde)                                   = res43;
//  local_residual.segment(_n_xi_Sorp_tilde,_n_xi_Min_tilde)                = res44;
//  local_residual.segment(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp)  = res45;
//  local_residual.tail(_n_xi_Min)                                          = res46;
//
//  for (i=0; i < _n_xi_global; i++)
//      _residual_global[i]->setValue(node_idx, local_residual[i]);
//
//
//  }
}


// element based assembly
template <class T_DIS_SYS, class T_USER_FUNCTION_DATA>
void TemplateTransientResidualFEMFunction_GIA_Reduct<T_DIS_SYS, T_USER_FUNCTION_DATA>::assembly(const NumLib::TimeStep & delta_t, const std::size_t _n_xi, const SolutionLib::SolutionVector & u_cur_xiglob,
                                            SolutionLib::SolutionVector & residual_global, std::size_t switchOn, MathLib::LocalMatrix &_mat_Asorp, MathLib::LocalMatrix &_mat_Amin)
        {
    //MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
    MeshLib::IMesh* msh = _dis_sys->getMesh();
    const size_t n_ele = msh->getNumberOfElements();
    double _theta(1.0);

    MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde, loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;

    // previous xi global
    MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;

    MathLib::LocalVector loc_vec_Rate, vec_Rate_temp;

    // initialize the local vector
    //current xi global
    loc_cur_xi_global               = MathLib::LocalVector::Zero( _n_xi_global );
    loc_cur_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde );
    loc_cur_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde );
    loc_cur_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp);
    loc_cur_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min);
    loc_cur_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin);
    //previous xi global
    loc_pre_xi_global               = MathLib::LocalVector::Zero( _n_xi_global );
    loc_pre_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde );
    loc_pre_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde );
    loc_pre_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp);
    loc_pre_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min);
    loc_pre_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin);


    for (size_t i=0; i<n_ele; i++)
    {
        MeshLib::IElement *e = msh->getElement(i);

        std::vector<size_t> ele_node_ids;
        //MathLib::Vector<size_t> ele_node_ids;
        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);

        //TODO get the water content and multiply it in LHS and RHS

        double dt = delta_t.getTimeStepSize();

        MathLib::LocalMatrix localM = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
        MathLib::LocalMatrix localK = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
        MathLib::LocalMatrix localDispersion = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
        MathLib::LocalMatrix localAdvection = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
        MathLib::LocalVector F = MathLib::LocalVector::Zero(ele_node_ids.size());
        MathLib::LocalMatrix dispersion_diffusion;
        MathLib::LocalMatrix d_poro = MathLib::LocalMatrix::Zero(3,3);
        MathLib::LocalMatrix poro(1,1);
        NumLib::ITXFunction::DataType v;

        //assembleODE(time, e, local_u_n1, local_u_n, M, K, F);
        _fe = _function_data->get_feObjects()->getFeObject(*e);

        const size_t n_dim = e->getDimension();
        size_t mat_id = e->getGroupID();
        MaterialLib::PorousMedia* _pm = Ogs6FemData::getInstance()->list_pm[mat_id];

        localDispersion.setZero(localK.rows(), localK.cols());
        localAdvection.setZero (localK.rows(), localK.cols());

        double cmp_mol_diffusion = .0;

        _q = _fe->getIntegrationMethod();
        double gp_x[3], real_x[3];
        poro.setZero();
        d_poro.setZero();
        double disp_l = 0.0;
        double disp_t = 0.0;

        for (size_t j=0; j < _q->getNumberOfSamplingPoints(); j++)
        {
            _q->getSamplingPoint(j, gp_x);
            _fe->computeBasisFunctions(gp_x);
            _fe->getRealCoordinates(real_x);
            NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e->getID(), j, real_x);

            _pm->porosity->eval(gp_pos, poro);
            _pm->dispersivity_long->eval(gp_pos, disp_l);
            _pm->dispersivity_trans->eval(gp_pos, disp_t);

            d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
            d_poro(1,1) = cmp_mol_diffusion * poro(0,0);
            d_poro(2,2) = cmp_mol_diffusion * poro(0,0);
            _function_data->_vel->eval(gp_pos, v);
            NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();

            // calculating dispersion tensor according to Benchmark book p219, Eq. 10.15
            // D_{ij} = \alpha_T |v| \delta_{ij} + (\alpha_L - \alpha_T) \frac{v_i v_j}{|v|} + D^{d}_{ii}
            dispersion_diffusion.setIdentity(n_dim, n_dim);
            dispersion_diffusion *= disp_l * v.norm();
            dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm();
            dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);
            // --------debugging--------------
            // std::cout << "dispersion_diffusion Matrix" << std::endl;
            // std::cout << dispersion_diffusion << std::endl;
            // --------end of debugging-------

            _fe->integrateWxN(j, poro, localM);
            _fe->integrateDWxDN(j, dispersion_diffusion, localDispersion);
            _fe->integrateWxDN(j, v2, localAdvection);
        }

        localK = localDispersion + localAdvection;

        // mass lumping----------------------------
        for (size_t idx_ml=0; idx_ml < localM.rows(); idx_ml++ )
        {
            double mass_lump_val;
            mass_lump_val = localM.row(idx_ml).sum();
            localM.row(idx_ml).setZero();
            localM(idx_ml, idx_ml) = mass_lump_val;
        }

        //std::cout << "M="; M.write(std::cout); std::cout << std::endl;
        //std::cout << "K="; K.write(std::cout); std::cout << std::endl;

        //MathLib::LocalVector localLHS = MathLib::LocalVector::Zero(ele_node_ids.size());
        //MathLib::LocalVector localRHS = MathLib::LocalVector::Zero(ele_node_ids.size());
        MathLib::LocalVector local_u1 = MathLib::LocalVector::Zero(ele_node_ids.size());
        MathLib::LocalVector local_u0 = MathLib::LocalVector::Zero(ele_node_ids.size());
        MathLib::LocalVector local_u1_tilde = MathLib::LocalVector::Zero(ele_node_ids.size());
        MathLib::LocalVector local_u0_tilde = MathLib::LocalVector::Zero(ele_node_ids.size());

        //MathLib::LocalVector node_indx = MathLib::LocalVector::Zero(_n_xi);
        std::size_t xi_count, node_indx(0);
        double localLHS, localRHS;

        for (xi_count = 0; xi_count < _n_xi; xi_count++)
        {

    //      for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
    //      {
    //
    //          node_indx           = ele_node_ids[idx] * _n_xi + xi_count;
    //          local_u1(idx)       = _global_u1_tilde[node_indx];
    //          local_u0(idx)       = _global_u0[node_indx];
    //          local_u1_tilde(idx) = _global_u1_tilde[node_indx];
    //          local_u0_tilde(idx) = _global_u0_tilde[node_indx];
    //      }


            for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
            {
                node_indx = ele_node_ids[idx] * _n_xi + xi_count;

                for (i=0; i < _n_xi_global; i++)
                    loc_cur_xi_global[i] = u_cur_xiglob[ele_node_ids[idx] * _n_xi_global + xi_count];
                for (i=0; i < _n_xi_global; i++)
                    loc_pre_xi_global[i] = this->_xi_global[i]->getValue(ele_node_ids[idx]);

                loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde);
                loc_cur_xi_Min_tilde    = loc_cur_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
                loc_cur_xi_Sorp         = loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
                loc_cur_xi_Min          = loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);

                loc_pre_xi_Sorp_tilde   = loc_pre_xi_global.head(_n_xi_Sorp_tilde);
                loc_pre_xi_Min_tilde    = loc_pre_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
                loc_pre_xi_Sorp         = loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
                loc_pre_xi_Min          = loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);

                for (i=0; i < _J_tot_kin; i++)
                    loc_vec_Rate[i] = this->_global_vec_Rate[i]->getValue(ele_node_ids[idx]);

                MathLib::LocalVector vec_Rate_temp;
                if(switchOn)  //if switchOn calculate for xi sorp tilde and for xi min tilde otherwise
                {
                                local_u1(idx)       = loc_cur_xi_Sorp[node_indx];
                                local_u0(idx)       = loc_pre_xi_Sorp[node_indx];
                                local_u1_tilde(idx) = loc_cur_xi_Sorp_tilde[node_indx];
                                local_u0_tilde(idx) = loc_pre_xi_Sorp_tilde[node_indx];
                                vec_Rate_temp = _mat_Asorp * loc_vec_Rate;

                }
                else
                {
                                local_u1(idx)       = loc_cur_xi_Min[node_indx];
                                local_u0(idx)       = loc_pre_xi_Min[node_indx];
                                local_u1_tilde(idx) = loc_cur_xi_Min_tilde[node_indx];
                                local_u0_tilde(idx) = loc_pre_xi_Min_tilde[node_indx];
                                vec_Rate_temp = _mat_Amin * loc_vec_Rate;
                }

            // LHS = dt*(1/dt M + theta K)
            //localLHS = (localM * global_u1_tilde(node_indx)) + (delta_t * _theta * localK * global_u1(node_indx));
             localLHS = (localM.row(idx).dot(local_u1_tilde)) + (dt * _theta * (localK.row(idx).dot(local_u1)));
            //_global_vec_LHS(node_indx) = _global_vec_LHS(node_indx) + localLHS;
             // RHS = (1/dt M - (1-theta) K) u0 + F
             localRHS = (localM.row(idx).dot(local_u0_tilde)) - (dt * (1.-_theta) * localK.row(idx).dot(local_u0)); //local_u_n;
             //_global_vec_RHS(node_indx) = _global_vec_RHS(node_indx) + localRHS;

            residual_global[_n_xi_global * ele_node_ids[idx] + xi_count] = localLHS - localRHS - vec_Rate_temp[xi_count];

            }
        }
    }
}
