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
          _xi_global_pre(function_data->get_xi_global_pre()), _xi_local_new(function_data->get_xi_local_new()), _eta(function_data->get_eta()), _eta_bar(function_data->get_eta_bar()),
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
	void GlobalResidualAssembler(const NumLib::TimeStep & delta_t, 
                                 const SolutionLib::SolutionVector & u_cur_xiglob, 
                                 const SolutionLib::SolutionVector & u_pre_xiglob, 
                                 SolutionLib::SolutionVector & residual_global);

	void assembly(const NumLib::TimeStep & delta_t, 
                  const SolutionLib::SolutionVector & u_cur_xiglob,
                  const SolutionLib::SolutionVector & u_pre_xiglob,
                  SolutionLib::SolutionVector & residual_global, 
                  MathLib::LocalMatrix &_mat_Asorp, 
                  MathLib::LocalMatrix &_mat_Amin);

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
    std::vector<MyNodalFunctionScalar*> &_xi_global_pre, &_xi_local_new, &_eta, &_eta_bar, &_global_vec_Rate;
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

    // assembly
    MeshLib::IMesh* msh = _dis_sys->getMesh();
    r = .0;
    //node based operations
	// calculate the global residual
    GlobalResidualAssembler(*this->_time_step, u_n1, *_u_n0, r);

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
void TemplateTransientResidualFEMFunction_GIA_Reduct<T_DIS_SYS, T_USER_FUNCTION_DATA>
     ::GlobalResidualAssembler(const NumLib::TimeStep & delta_t,
                               const SolutionLib::SolutionVector & u_cur_xiglob, 
                               const SolutionLib::SolutionVector & u_pre_xiglob, 
                               SolutionLib::SolutionVector & residual_global)
{
    const size_t nnodes = _dis_sys->getMesh()->getNumberOfNodes();
    const double theta_water_content = 0.5;
    size_t j; 
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
    MathLib::LocalVector res43, res44, res45, res46, res47, local_residual;//, residual_global_res43, residual_global_res44, residual_global_res45, residual_global_res46, global_vec_LHS_sorp, global_vec_RHS_sorp, global_vec_LHS_min, global_vec_RHS_min;

    // eta vectors
    MathLib::LocalVector loc_cur_eta, loc_cur_eta_bar;

    // rate vector
    MathLib::LocalVector vec_Rate, vec_Rate_45, vec_Rate_46; //, global_vec_Rate_45, global_vec_Rate_46;

    MathLib::LocalMatrix mat_A1sorp, mat_A2sorpli, mat_A2sorpld, mat_A1min, mat_Ald, mat_A1kin;

    // initialize the local vector
    //current xi global
    loc_cur_xi_global               = MathLib::LocalVector::Zero( _n_xi_global );
    loc_cur_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde );
    loc_cur_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde );
    loc_cur_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp);
    loc_cur_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min);
    loc_cur_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin);

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
    vec_Rate					= MathLib::LocalVector::Zero( _J_tot_kin);

    // current eta mobie and immobile
    loc_cur_eta             = MathLib::LocalVector::Zero( _n_eta );
    loc_cur_eta_bar         = MathLib::LocalVector::Zero( _n_eta_bar );

    //initialize it before?
    mat_A1sorp        = _ReductionGIA->get_matrix_A1sorp();
    mat_A2sorpli      = _ReductionGIA->get_matrix_A2sorpli();
    mat_A1min         = _ReductionGIA->get_matrix_A1min();
    mat_Ald           = _ReductionGIA->get_matrix_Ald();
    mat_A2sorpld      = _ReductionGIA->get_matrix_A2sorpld();
    mat_A1kin         = _ReductionGIA->get_matrix_A1kin(); 

    MathLib::LocalMatrix mat_Asorp = mat_A1sorp - mat_A2sorpli;
    MathLib::LocalMatrix mat_Amin  = mat_A1min  - mat_Ald * mat_A2sorpld;

    // loop over all the nodes
    for (size_t node_idx = 0; node_idx < nnodes; node_idx++ )
    {

            // on each node, get the right start value
            // get the right set of xis

            for (size_t i=0; i < _n_xi_global; i++)
                loc_cur_xi_global[i] = u_cur_xiglob[node_idx * _n_xi_global + i];
            for (size_t i=0; i < _n_xi_local; i++)
                loc_cur_xi_local[i] = _xi_local_new[i]->getValue(node_idx);
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

            if(_n_xi_Kin != 0){
            // calculate the nodal kinetic reaction rates
            _ReductionGIA->Calc_Kin_Rate_temp(loc_cur_xi_Mob,
                                         loc_cur_xi_Sorp,
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

            res45 = delta_t.getTimeStepSize() * theta_water_content * mat_Asorp * vec_Rate; 
            res46 = delta_t.getTimeStepSize() * theta_water_content * mat_Amin  * vec_Rate; 
            res47 = delta_t.getTimeStepSize() * theta_water_content * mat_A1kin * vec_Rate; 

            for (j=0; j<_n_xi_Sorp_tilde; j++ )
                residual_global[_n_xi_global * node_idx + _n_xi_Sorp_tilde + _n_xi_Min_tilde + j] -= res45(j) ;
            for (j=0; j<_n_xi_Min_tilde; j++ )
                residual_global[_n_xi_global * node_idx + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp_tilde + j] -= res46(j) ;
            for (j=0; j<_n_xi_Kin; j++ )
                residual_global[_n_xi_global * node_idx + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp_tilde + _n_xi_Min_tilde + j] -= res47(j) ;
            }
     
    } // end of loop over all nodes. 

    // loop over all elements and assemble Eq. 45-47
    this->assembly(delta_t, u_cur_xiglob, u_pre_xiglob, residual_global, mat_Asorp, mat_Amin);

}


// element based assembly
template <class T_DIS_SYS, class T_USER_FUNCTION_DATA>
void TemplateTransientResidualFEMFunction_GIA_Reduct
     <T_DIS_SYS, T_USER_FUNCTION_DATA>
     ::assembly(const NumLib::TimeStep & delta_t, 
                const SolutionLib::SolutionVector & u_cur_xiglob,
                const SolutionLib::SolutionVector & u_pre_xiglob,
                SolutionLib::SolutionVector & residual_global,  
                MathLib::LocalMatrix &_mat_Asorp, 
                MathLib::LocalMatrix &_mat_Amin)
{
    //MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
    MeshLib::IMesh* msh = _dis_sys->getMesh();
    const size_t n_ele = msh->getNumberOfElements();
    const size_t n_max_connect_nodes = 20; 
    size_t nnodes; 
    double _theta(1.0);
    size_t _n_xi_trans = this->_n_xi_Sorp + this->_n_xi_Min + this->_n_xi_Kin; 
    size_t i,j,k; 

    MathLib::LocalVector loc_cur_xi_global, 
                         loc_cur_xi_Sorp_tilde, 
                         loc_cur_xi_Min_tilde, 
                         loc_cur_xi_Sorp, 
                         loc_cur_xi_Min, 
                         loc_cur_xi_Kin;

    // previous xi global
    MathLib::LocalVector loc_pre_xi_global, 
                         loc_pre_xi_Sorp_tilde, 
                         loc_pre_xi_Min_tilde, 
                         loc_pre_xi_Sorp, 
                         loc_pre_xi_Min, 
                         loc_pre_xi_Kin;

    MathLib::LocalVector loc_vec_Rate, vec_Rate_temp;

    MathLib::LocalVector localLHS_xi_sorp,
                         localRHS_xi_sorp, 
                         localLHS_xi_min, 
                         localRHS_xi_min,
                         localLHS_xi_kin, 
                         localRHS_xi_kin, 
                         local_res_sorp, 
                         local_res_min,
                         local_res_kin; 

    loc_vec_Rate  = MathLib::LocalVector::Zero( _J_tot_kin );
    vec_Rate_temp = MathLib::LocalVector::Zero( _J_tot_kin );

    // initialize the local vector
    //current xi global
    loc_cur_xi_global               = MathLib::LocalVector::Zero( _n_xi_global * n_max_connect_nodes );
    loc_cur_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde * n_max_connect_nodes  );
    loc_cur_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde * n_max_connect_nodes  );
    loc_cur_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp * n_max_connect_nodes );
    loc_cur_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min * n_max_connect_nodes );
    loc_cur_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin * n_max_connect_nodes );
    //previous xi global
    loc_pre_xi_global               = MathLib::LocalVector::Zero( _n_xi_global * n_max_connect_nodes  );
    loc_pre_xi_Sorp_tilde           = MathLib::LocalVector::Zero( _n_xi_Sorp_tilde * n_max_connect_nodes  );
    loc_pre_xi_Min_tilde            = MathLib::LocalVector::Zero( _n_xi_Min_tilde * n_max_connect_nodes  );
    loc_pre_xi_Sorp                 = MathLib::LocalVector::Zero( _n_xi_Sorp * n_max_connect_nodes );
    loc_pre_xi_Min                  = MathLib::LocalVector::Zero( _n_xi_Min * n_max_connect_nodes );
    loc_pre_xi_Kin                  = MathLib::LocalVector::Zero( _n_xi_Kin * n_max_connect_nodes );


    for ( i=0; i<n_ele; i++)
    {
        MeshLib::IElement *e = msh->getElement(i);

        std::vector<size_t> ele_node_ids;
        //MathLib::Vector<size_t> ele_node_ids;
        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
        // number of connecting nodes
        nnodes = ele_node_ids.size(); 

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

        //Local RHS, LHS and residual
        localLHS_xi_sorp  = MathLib::LocalVector::Zero(ele_node_ids.size());
        localRHS_xi_sorp  = MathLib::LocalVector::Zero(ele_node_ids.size());
        localLHS_xi_min   = MathLib::LocalVector::Zero(ele_node_ids.size());
        localRHS_xi_min   = MathLib::LocalVector::Zero(ele_node_ids.size());
        localLHS_xi_kin   = MathLib::LocalVector::Zero(ele_node_ids.size());
        localRHS_xi_kin   = MathLib::LocalVector::Zero(ele_node_ids.size());
        local_res_sorp    = MathLib::LocalVector::Zero(ele_node_ids.size());
        local_res_min     = MathLib::LocalVector::Zero(ele_node_ids.size());
        local_res_kin     = MathLib::LocalVector::Zero(ele_node_ids.size());


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

        for (j=0; j < _q->getNumberOfSamplingPoints(); j++)
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
            _vel->eval(gp_pos, v);
            NumLib::ITXFunction::DataType v2 = v.topRows(n_dim).transpose();

            // calculating dispersion tensor according to Benchmark book p219, Eq. 10.15
            // D_{ij} = \alpha_T |v| \delta_{ij} + (\alpha_L - \alpha_T) \frac{v_i v_j}{|v|} + D^{d}_{ii}
            dispersion_diffusion.setIdentity(n_dim, n_dim);

            dispersion_diffusion *= disp_l * v.norm();

            dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm();

            dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);

            _fe->integrateWxN(j, poro, localM);
            _fe->integrateDWxDN(j, dispersion_diffusion, localDispersion);
            _fe->integrateWxDN(j, v2, localAdvection);
        }

        //// --------debugging--------------
        //std::cout << "localDispersion Matrix" << std::endl;
        //std::cout << localDispersion << std::endl;
        //std::cout << "localAdvection Matrix" << std::endl;
        //std::cout << localAdvection << std::endl;
        //// --------end of debugging-------

        localK = localDispersion + localAdvection;

        // mass lumping----------------------------
        for (int idx_ml=0; idx_ml < localM.rows(); idx_ml++ )
        {
            double mass_lump_val;
            mass_lump_val = localM.row(idx_ml).sum();
            localM.row(idx_ml).setZero();
            localM(idx_ml, idx_ml) = mass_lump_val;
        }

        // --------debugging--------------
        //std::cout << "localK Matrix" << std::endl;
        //std::cout << localK << std::endl;
        //std::cout << "localM Matrix" << std::endl;
        //std::cout << localM << std::endl;
        // --------end of debugging-------

        //MathLib::LocalVector node_indx = MathLib::LocalVector::Zero(_n_xi);
        std::size_t idx_xi, node_indx, val_idx;

        for (idx_xi = 0; idx_xi < _n_xi_global ; idx_xi++)
        {
            for( node_indx = 0; node_indx < nnodes; node_indx++)
            {
                // calculating the global index 
                val_idx = ele_node_ids[node_indx] * _n_xi_global + idx_xi;
                // now we get the nodal values connected to this element
                loc_cur_xi_global[idx_xi*nnodes+node_indx] = u_cur_xiglob[val_idx];
                loc_pre_xi_global[idx_xi*nnodes+node_indx] = u_pre_xiglob[val_idx];
            }  // end of idx_xi
        }  // end of node_idx

        // --------debugging--------------
//        std::cout << "loc_cur_xi_global" << std::endl;
//        std::cout << loc_cur_xi_global << std::endl;
//        std::cout << "loc_pre_xi_global" << std::endl;
//        std::cout << loc_pre_xi_global << std::endl;
        // --------end of debugging-------

        loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde*nnodes);
        loc_cur_xi_Min_tilde    = loc_cur_xi_global.segment(_n_xi_Sorp_tilde*nnodes, _n_xi_Min_tilde*nnodes);
        loc_cur_xi_Sorp         = loc_cur_xi_global.segment(_n_xi_Sorp_tilde*nnodes  + _n_xi_Min_tilde*nnodes, _n_xi_Sorp*nnodes);
        loc_cur_xi_Min          = loc_cur_xi_global.segment( _n_xi_Sorp_tilde*nnodes + _n_xi_Min_tilde*nnodes + _n_xi_Sorp*nnodes, _n_xi_Min*nnodes);
        loc_cur_xi_Kin          = loc_cur_xi_global.segment( _n_xi_Sorp_tilde*nnodes + _n_xi_Min_tilde*nnodes + _n_xi_Sorp*nnodes + _n_xi_Min*nnodes, _n_xi_Kin*nnodes);

        loc_pre_xi_Sorp_tilde   = loc_pre_xi_global.head(_n_xi_Sorp_tilde*nnodes);
        loc_pre_xi_Min_tilde    = loc_pre_xi_global.segment(_n_xi_Sorp_tilde*nnodes, _n_xi_Min_tilde*nnodes);
        loc_pre_xi_Sorp         = loc_pre_xi_global.segment(_n_xi_Sorp_tilde*nnodes  + _n_xi_Min_tilde*nnodes, _n_xi_Sorp*nnodes);
        loc_pre_xi_Min          = loc_pre_xi_global.segment( _n_xi_Sorp_tilde*nnodes + _n_xi_Min_tilde*nnodes + _n_xi_Sorp*nnodes, _n_xi_Min*nnodes);
        loc_pre_xi_Kin          = loc_pre_xi_global.segment( _n_xi_Sorp_tilde*nnodes + _n_xi_Min_tilde*nnodes + _n_xi_Sorp*nnodes + _n_xi_Min*nnodes, _n_xi_Kin*nnodes);

        //    --------debugging--------------
//        std::cout << "loc_cur_xi_Kin" << std::endl;
//        std::cout << loc_cur_xi_Kin << std::endl;
//        std::cout << "loc_pre_xi_Kin" << std::endl;
//        std::cout << loc_pre_xi_Kin << std::endl;
        //    --------end of debugging-------

        // LHS = (1/dt M + theta theta K) u1
        // RHS = (1/dt M - (1-theta) K) u0 + F
        for ( j=0; j < _n_xi_Sorp_tilde; j++ )
        {
            localLHS_xi_sorp  = localM * loc_cur_xi_Sorp_tilde.segment(j*nnodes, nnodes) 
                                  + dt * _theta * localK * loc_cur_xi_Sorp.segment(j*nnodes, nnodes);

            localRHS_xi_sorp  = localM * loc_pre_xi_Sorp_tilde.segment(j*nnodes, nnodes) 
                                  + dt * (1.0 - _theta) * localK * loc_pre_xi_Sorp.segment(j*nnodes, nnodes);
            local_res_sorp    = localLHS_xi_sorp - localRHS_xi_sorp;
            for (k=0; k<nnodes; k++)
            {
                residual_global[_n_xi_global * ele_node_ids[k] + _n_xi_Sorp_tilde + _n_xi_Min_tilde + j] += local_res_sorp(k) ;
            }
        }  // end of for j
        for ( j=0; j <_n_xi_Min_tilde; j++ )
        {
            localLHS_xi_min  = localM * loc_cur_xi_Min_tilde.segment(j*nnodes, nnodes) 
                                 + dt * _theta * localK * loc_cur_xi_Min.segment(j*nnodes, nnodes);

            // --------debugging--------------
            //MathLib::LocalVector tempp_tilde = loc_cur_xi_Min_tilde.segment(j*nnodes, nnodes);
            //MathLib::LocalVector tempp = loc_cur_xi_Min_tilde.segment(j*nnodes, nnodes);
            //std::cout << "loc_cur_xi_Min_tilde" << std::endl;
            //std::cout << tempp_tilde << std::endl;
            //
            //std::cout << "loc_cur_xi_Min" << std::endl;
            //std::cout << tempp << std::endl;
            //std::cout << "loc_cur_xi_Min_tilde" << std::endl;
            //std::cout << tempp << std::endl;
            // --------end of debugging-------

            localRHS_xi_min  = localM * loc_pre_xi_Min_tilde.segment(j*nnodes, nnodes) 
                                 + dt * (1.0 - _theta) * localK * loc_pre_xi_Min.segment(j*nnodes, nnodes);
            local_res_min    = localLHS_xi_min - localRHS_xi_min; 
            // --------debugging--------------
            //std::cout << "localLHS_xi_min vector" << std::endl;
            //std::cout << localLHS_xi_min << std::endl;
            //std::cout << "localRHS_xi_min vector" << std::endl;
            //std::cout << localRHS_xi_min << std::endl;
            //std::cout << "local_res_min vector" << std::endl;
            //std::cout << local_res_min << std::endl;
            // --------end of debugging-------
            for (k=0; k<nnodes; k++)
            {
                residual_global[_n_xi_global * ele_node_ids[k] + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + j] += local_res_min(k) ;
                // --------debugging--------------
                //double tmp = _n_xi_global * ele_node_ids[k] + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + j;
                //std::cout << "global nodal indx:" << std::endl;
                //std::cout << tmp << std::endl;
                // --------end of debugging-------
            }
        }  // end of for j
        for ( j=0; j <_n_xi_Kin; j++ )
        {
            localLHS_xi_kin  = localM * loc_cur_xi_Kin.segment(j*nnodes, nnodes) 
                                 + dt * _theta * localK * loc_cur_xi_Kin.segment(j*nnodes, nnodes);
            localRHS_xi_kin  = localM * loc_pre_xi_Kin.segment(j*nnodes, nnodes) 
                                 + dt * (1.0 - _theta) * localK * loc_pre_xi_Kin.segment(j*nnodes, nnodes);
            local_res_kin    = localLHS_xi_kin - localRHS_xi_kin; 
            for (k=0; k<nnodes; k++)
            {
                residual_global[_n_xi_global * ele_node_ids[k] + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min + j] += local_res_kin(k) ;
            }
        }  // end of for j
    }  // end of loop over all elements
}
