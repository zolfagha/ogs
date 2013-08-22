/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemDxEQS.h
 *
 * Created on 2013-08-016 by Reza Zolfaghari & Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "MathLib/DataType.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TransientAssembler/TransientElementWiseMatrixUpdater.h"
#include "FemLib/Function/FemFunction.h"
#include "SolutionLib/DataType.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"
#include "SolutionLib/Fem/FemVariable.h"
#include "IterativeLinearSolvers"

/**
 * \brief Template class for transient linear FEM functions
 *
 * \tparam T_SPACE_ASSEMBLER
 */
template <
    class T_DIS_SYS,
    class T_LINEAR_SOLVER,
    class T_USER_FUNCTION_DATA
//    class T_LOCAL_JACOBIAN_ASSEMBLER
    >
class TemplateTransientDxFEMFunction_GIA_Reduct
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef typename T_DIS_SYS::template MyLinearEquation<T_LINEAR_SOLVER,DiscreteLib::SparsityBuilderFromNodeConnectivity>::type MyLinaerEQS;
    typedef SolutionLib::FemVariable MyVariable;
    typedef T_LINEAR_SOLVER LinearSolverType;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
//    typedef T_LOCAL_JACOBIAN_ASSEMBLER UserLocalJacobianAssembler;
//    typedef typename NumLib::TransientElementWiseMatrixUpdater<UserLocalJacobianAssembler> MyUpdater;
//    typedef typename T_DIS_SYS::template MyLinearEquationAssembler<MyUpdater,LinearSolverType>::type MyGlobalAssembler;

    /// constructor
    /// @param problem        Fem problem
    /// @param linear_eqs    Discrete linear equation
    TemplateTransientDxFEMFunction_GIA_Reduct(MeshLib::IMesh* msh,
    		                                  std::vector<MyVariable*> &list_var,
    		                                  MyLinaerEQS* linear_eqs,
    		                                  DiscreteLib::DofEquationIdTable* dofManager,
    		                                  T_USER_FUNCTION_DATA* userData)
        : _msh(msh),  _linear_eqs(linear_eqs),
          _t_n1(0), _u_n0(0), _list_var(list_var), _dofManager(dofManager), _userData(userData), _ReductionGIA(userData->getReductionGIA()),
     	 _n_Comp(_ReductionGIA->get_n_Comp()), _I_mob(_ReductionGIA->get_n_Comp_mob()), _I_min(_ReductionGIA->get_n_Comp_min()), _n_xi_Kin_bar(_ReductionGIA->get_n_xi_Kin_bar())
     	, _n_xi_Mob(_ReductionGIA->get_n_xi_Mob()), _n_eta(_ReductionGIA->get_n_eta()), _n_eta_bar(_ReductionGIA->get_n_eta_bar()), _n_xi_Sorp_tilde(_ReductionGIA->get_n_xi_Sorp_tilde()),
     	_n_xi_Min_tilde(_ReductionGIA->get_n_xi_Min_tilde()), _n_xi_Sorp(_ReductionGIA->get_n_xi_Sorp()), _n_xi_Min(_ReductionGIA->get_n_xi_Min()), _n_xi_Sorp_bar_li(_ReductionGIA->get_n_xi_Sorp_bar_li()),
     	_n_xi_Sorp_bar_ld(_ReductionGIA->get_n_xi_Sorp_bar_ld()), _n_xi_Kin(_ReductionGIA->get_n_xi_Kin()), _n_xi_Min_bar(_ReductionGIA->get_n_xi_Min_bar()), _I_NMin_bar(_ReductionGIA->get_n_Comp_NMin_bar())
    {
    };

    ///
    virtual ~TemplateTransientDxFEMFunction_GIA_Reduct() {};

    void setVelocity(const NumLib::ITXFunction *vel)
    {
        _vel = const_cast<NumLib::ITXFunction*>(vel);
    }

    ///
    NumLib::TemplateFunction<SolutionLib::SolutionVector,SolutionLib::SolutionVector>* clone() const
    {
        return new TemplateTransientDxFEMFunction_GIA_Reduct<
                    T_DIS_SYS,
                    T_LINEAR_SOLVER,
                    T_USER_FUNCTION_DATA
                    >(_list_var, _linear_eqs, _dofManager, _userData);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const SolutionLib::SolutionVector* u_n0)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
    };

    /// solve Newton-Raphson equations
    /// @param u_n1    current solution
    /// @param r     residual
    /// @param du    increment
    void eval(const SolutionLib::SolutionVector &u_n1, const SolutionLib::SolutionVector &r, SolutionLib::SolutionVector &du);

private:
    void GlobalJacobianAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, LinearSolverType & eqsJacobian_global);
    void Vprime( MathLib::LocalVector & vec_conc, MathLib::LocalVector & logk_min, MathLib::LocalMatrix & mat_S1min, MathLib::LocalMatrix & mat_S1mob, MathLib::LocalMatrix & mat_S1sorp, MathLib::LocalMatrix & mat_S1sorpli,
    		 MathLib::LocalMatrix & mat_S1kin_ast, MathLib::LocalMatrix & mat_S2sorp, MathLib::LocalMatrix & mat_vprime);
    void NumDiff(std::size_t & col, const double & delta_xi, ogsChem::LocalVector & f, ogsChem::LocalVector & f_old, ogsChem::LocalVector & unknown, ogsChem::LocalVector & DrateDxi);
    void AddMassLaplasTerms(const double & dt, const double & theta_water_content, LinearSolverType & eqsJacobian_global);
    void cal_nodal_rate(ogsChem::LocalVector &xi, ogsChem::LocalVector &local_eta_bar, ogsChem::LocalVector &local_eta, std::size_t & n_xi_Kin_total, ogsChem::LocalMatrix & mat_S1_ast, ogsChem::LocalMatrix & mat_S2_ast,
    					ogsChem::LocalMatrix & mat_S1_orth, ogsChem::LocalMatrix & mat_S2_orth, ogsChem::LocalVector &vec_rate_new);
private:
    MeshLib::IMesh* _msh;
    DiscreteLib::DofEquationIdTable* _dofManager;
    MyLinaerEQS* _linear_eqs;
    const NumLib::TimeStep* _t_n1;
    const SolutionLib::SolutionVector* _u_n0;
    std::vector<MyVariable*> _list_var;
    T_USER_FUNCTION_DATA* _userData;

    //TODO pass via constructor
    ogsChem::chemReductionGIA* _ReductionGIA;
    std::map<size_t, ReductionGIANodeInfo*>* _bc_info;
    std::vector<MyNodalFunctionScalar*> _concentrations, _xi_global, _xi_local, _eta, _eta_bar, _global_vec_Rate;
    std::vector<ogsChem::chemReactionKin*>  _list_kin_reactions;
    NumLib::ITXFunction* _vel;
    //TODO set the followings from _ReductionGIA
    size_t _n_xi_global, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Sorp, _n_xi_Min, _n_xi_Kin, _n_xi_local, _n_xi_Sorp_bar, _n_xi_Min_bar, _n_eta, _n_eta_bar, _n_xi_Mob, _n_xi_Kin_bar, _vec_Rate_rows;
    size_t _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_ld, _n_Comp, _I_NMin_bar, _I_mob, _I_min;
};


template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::eval(const SolutionLib::SolutionVector &u_n1, const SolutionLib::SolutionVector &r, SolutionLib::SolutionVector &du)
{
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
    const SolutionLib::SolutionVector* u_n = this->_u_n0;

    _linear_eqs->initialize();

    // setup BC1 for solution increment
    for (size_t i=0; i<_list_var.size(); i++) {
        MyVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
            SolutionLib::FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            std::vector<double> bc_value_for_dx(bc1->getListOfBCNodes().size(), .0);
            for (size_t k=0; k<bc_value_for_dx.size(); k++) {
                size_t id = bc1->getListOfBCNodes()[k];
                double v =  bc1->getListOfBCValues()[k];
                size_t eqs_id = _linear_eqs->getDofMapManger()->mapEqsID(i, _msh->getID(), id);
                bc_value_for_dx[k] = v - u_n1[eqs_id]; // dx = (bc value) - (current value)
            }
            var_bc_id.insert(var_bc_id.end(), bc1->getListOfBCNodes().begin(), bc1->getListOfBCNodes().end());
            var_bc_val.insert(var_bc_val.end(), bc_value_for_dx.begin(), bc_value_for_dx.end());
        }
        _linear_eqs->setPrescribedDoF(i, var_bc_id, var_bc_val);
    }

    // TODO global assembly
    //NumLib::ElementWiseTransientDxEQSAssembler<UserLocalJacobianAssembler> assembler(&t_n1, u_n, &u_n1, _local_assembler);
//    MyUpdater updater(&t_n1, _msh, u_n, &u_n1, _local_assembler);
//    MyGlobalAssembler assembler(&updater);
//    _linear_eqs->construct(assembler);
    LinearSolverType* eqsJacobian_global = _linear_eqs->getLinearSolver();
    GlobalJacobianAssembler(t_n1.getTimeStepSize(), u_n1, *eqsJacobian_global);

    // set residual
    _linear_eqs->addRHS(r, -1.0);

    // solve
    _linear_eqs->solve();
    _linear_eqs->getX(du);
}

template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::GlobalJacobianAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, LinearSolverType & eqsJacobian_global)
{
    //---------------------------------------------------------------------------
    // Note for Reza:
    // eqsJacobian_global is actually MathLib/LinAlg/LisLinearEquation object in your case
    // In this function, you have to update a coefficient matrix in the linear equation object.
    // For example, you can add a value to the jacobian matrix as
    //
    // eqsJacobian_global.addA(0,0, 1.0);
    //
    // You can get a dimension of the matrix as
    //
    // size_t eqs_dim = eqsJacobian_global.getDimension();
    //
    //---------------------------------------------------------------------------


    using namespace std::placeholders;
    size_t i, node_idx, indx_tmp, nnodes;
    _n_xi_Sorp_bar_li  = _ReductionGIA->get_n_xi_Sorp_bar_li();
    _n_xi_Sorp_bar_ld  = _ReductionGIA->get_n_xi_Sorp_bar_ld();
    MyDiscreteSystem* dis = 0;
    nnodes = dis->getMesh()->getNumberOfNodes();
    const double theta_water_content(0.32);
    const double delta_xi = 1E-14;

    // current xi global
    MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde,
                         loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;
    // previous xi global
    MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;

    // current xi local
    MathLib::LocalVector loc_cur_xi_local, loc_xi_local, loc_cur_xi_Mob, loc_cur_xi_Sorp_bar, loc_cur_xi_Min_bar, loc_cur_xi_Kin_bar, loc_cur_xi_Sorp_bar_li, loc_cur_xi_Sorp_bar_ld;

    // eta vectors
    MathLib::LocalVector loc_cur_eta, loc_cur_eta_bar;

    // local concentration vector
    MathLib::LocalVector vec_conc;

    MathLib::LocalMatrix mat_p1Fder,mat_p1Ftrans, mat_p1F, mat_p2F, mat_Global_Jacobian;
    std::size_t n_xi_Kin_total = _ReductionGIA->get_n_xi_Kin_total();

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
    // current eta mobie and immobile
    loc_cur_eta             = MathLib::LocalVector::Zero( _n_eta );
    loc_cur_eta_bar         = MathLib::LocalVector::Zero( _n_eta_bar );
    vec_conc                = MathLib::LocalVector::Zero(_n_Comp);

    MathLib::LocalMatrix mat_vprime  = MathLib::LocalMatrix::Zero(_n_xi_local,_n_xi_global);
    MathLib::LocalMatrix Jacobian_local = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    MathLib::LocalVector vec_rate_old = MathLib::LocalVector::Zero(_vec_Rate_rows);
    MathLib::LocalVector vec_Rate = MathLib::LocalVector::Zero(_vec_Rate_rows);
    MathLib::LocalVector logk_min;
    MathLib::LocalMatrix mat_S1min, mat_S1mob, mat_S1sorp, mat_S1sorpli, mat_S1kin_ast, mat_S2sorp;

    mat_p1Fder  = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p1Ftrans = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p1F     = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p2F     = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    double dt = delta_t.getTimeStepSize();

    logk_min 		   =  _ReductionGIA->get_logk_min();
    mat_S1min 		   =  _ReductionGIA->get_matrix_S1min();
    mat_S1mob 		   =  _ReductionGIA->get_matrix_S1mob();
    mat_S1sorp 		   =  _ReductionGIA->get_matrix_S1sorp();
    mat_S1sorpli 	   =  _ReductionGIA->get_matrix_S1sorpli();
    mat_S1kin_ast 	   =  _ReductionGIA->get_matrix_S1kin_ast();
    mat_S2sorp 		   =  _ReductionGIA->get_matrix_S2sorp();

    MathLib::LocalMatrix mat_Asorp, mat_Amin, mat_A1sorp, mat_A2sorpli, mat_A1min, mat_Ald, mat_A2sorpld;
    mat_A1sorp		   = _ReductionGIA->get_matrix_A1sorp();
    mat_A2sorpli 	   = _ReductionGIA->get_matrix_A2sorpli();
    mat_A1min		   = _ReductionGIA->get_matrix_A1min();
    mat_Ald 		   = _ReductionGIA->get_matrix_Ald();
    mat_A2sorpld  	   = _ReductionGIA->get_matrix_A2sorpld();
    MathLib::LocalMatrix mat_A1kin   = _ReductionGIA->get_matrix_A1kin();

    mat_Asorp = mat_A1sorp - mat_A2sorpli;
    mat_Amin = mat_A1min - mat_Ald * mat_A2sorpld;

    MathLib::LocalMatrix mat_S1_ast  = _ReductionGIA->get_mat_S1_ast();
    MathLib::LocalMatrix mat_S2_ast  = _ReductionGIA->get_mat_S2_ast();
    MathLib::LocalMatrix mat_S1_orth  = _ReductionGIA->get_mat_S1_orth();
    MathLib::LocalMatrix mat_S2_orth  = _ReductionGIA->get_mat_S2_orth();

    std::size_t n_xi_total = _n_xi_local + _n_xi_global;
    // loop over all the nodes
    for (size_t node_idx=0; node_idx < nnodes; node_idx++ )
    {
        // skip the boundary nodes
        if (  _bc_info->find( node_idx ) == _bc_info->end() )
        {

            // on each node, get the right start value
            // get the right set of xis

            for (i=0; i < _n_xi_global; i++)
                loc_cur_xi_global[i] = u_cur_xiglob[node_idx * _n_xi_global + i];
            for (i=0; i < _n_xi_local; i++)
                loc_cur_xi_local[i] = this->_xi_local[i]->getValue(node_idx);
            for (i=0; i < _n_eta; i++)
                loc_cur_eta[i] = this->_eta[i]->getValue(node_idx);
            // fill in eta_immob
            for (i=0; i < _n_eta_bar; i++)
                loc_cur_eta_bar[i] = this->_eta_bar[i]->getValue(node_idx);

            for (i=0; i < _vec_Rate_rows; i++)
                vec_rate_old[i] = this->_global_vec_Rate[i]->getValue(node_idx);

            for (i=0; i < _n_Comp; i++)
                vec_conc[i] = this->_concentrations[i]->getValue(node_idx);

//            // current xi global
//            loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde);
//            loc_cur_xi_Min_tilde    = loc_cur_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
//            loc_cur_xi_Sorp         = loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
//            loc_cur_xi_Min          = loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);
//            loc_cur_xi_Kin          = loc_cur_xi_global.tail(_n_xi_Kin);
//
//            // current xi local
//            loc_cur_xi_Mob          = loc_cur_xi_local.head(_n_xi_Mob);
//            loc_cur_xi_Sorp_bar     = loc_cur_xi_local.segment(_n_xi_Mob, _n_xi_Sorp_bar);
//            loc_cur_xi_Min_bar      = loc_cur_xi_local.segment( _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar);
//            loc_cur_xi_Kin_bar      = loc_cur_xi_local.tail(_n_xi_Kin_bar);
//
//            loc_cur_xi_Sorp_bar_li  = loc_cur_xi_Sorp_bar.topRows(_n_xi_Sorp_bar_li);
//            loc_cur_xi_Sorp_bar_ld  = loc_cur_xi_Sorp_bar.bottomRows(_n_xi_Sorp_bar_ld);


            ///// calculate partial1F

            // calculate the derivative
            MathLib::LocalVector Unknown_vec 	 = MathLib::LocalVector::Zero(n_xi_total);
            MathLib::LocalVector DrateDxi    	 = MathLib::LocalVector::Zero(n_xi_total);
            MathLib::LocalVector vec_rate_new    = MathLib::LocalVector::Zero(_vec_Rate_rows);
            double temp_xi;
            Unknown_vec.head(_n_xi_global) = loc_cur_xi_global;
            Unknown_vec.tail(_n_xi_local) = loc_cur_xi_local;

            for(std::size_t i = 0; i < n_xi_total ; i++ ){
            	ogsChem::LocalVector xi = ogsChem::LocalVector::Zero(n_xi_total);
            	xi        = Unknown_vec;
            	temp_xi = xi(i);
            	xi(i)     = xi(i) + (delta_xi * temp_xi);
            	cal_nodal_rate(xi, loc_cur_eta_bar, loc_cur_eta, n_xi_Kin_total, mat_S1_ast, mat_S2_ast, mat_S1_orth ,mat_S2_orth, vec_rate_new);
            	//DrateDxi(i) = ( vec_rate_new - vec_rate_old) / (delta_xi * temp_xi);
            }

            ogsChem::LocalVector der_sorpT_R =  MathLib::LocalVector::Zero(_n_xi_Sorp_tilde);
            der_sorpT_R 					 =  DrateDxi.head(_n_xi_Sorp_tilde);
            ogsChem::LocalVector der_minT_R  =  MathLib::LocalVector::Zero(_n_xi_Min_tilde);
            der_minT_R 				    	 =  DrateDxi.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
            ogsChem::LocalVector der_kin_R   =  MathLib::LocalVector::Zero(_n_xi_Kin);
            der_kin_R 				    	 =  DrateDxi.segment(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin);
            ogsChem::LocalVector der_mob_R   = MathLib::LocalVector::Zero(_n_xi_Mob);
            der_mob_R 				    	 =  DrateDxi.segment(_n_xi_global, _n_xi_Mob);
            ogsChem::LocalVector der_sorpB_R = MathLib::LocalVector::Zero(_n_xi_Sorp_bar);
            der_sorpB_R 				     =  DrateDxi.segment(_n_xi_global + _n_xi_Mob, _n_xi_Sorp_bar);
            ogsChem::LocalVector der_minB_R  = MathLib::LocalVector::Zero(_n_xi_Min_bar);
            der_minB_R	 				     =  DrateDxi.segment(_n_xi_global + _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar);



            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Sorp, _n_xi_Sorp_tilde)                         =  mat_Asorp * der_sorpT_R;
            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, 0, _n_xi_Min, _n_xi_Sorp_tilde)             =  mat_Amin * der_sorpT_R;
            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, 0, _n_xi_Kin, _n_xi_Sorp_tilde) =  mat_A1kin * der_sorpT_R;

            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde, _n_xi_Sorp, _n_xi_Min_tilde)                         =  mat_Asorp * der_minT_R;
            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Sorp_tilde, _n_xi_Min, _n_xi_Min_tilde)             =  mat_Amin * der_minT_R;
            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp_tilde, _n_xi_Kin, _n_xi_Min_tilde) =  mat_A1kin * der_minT_R;

            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp, _n_xi_Kin)                         =  mat_Asorp * der_kin_R;
            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Min, _n_xi_Kin)             =  mat_Amin * der_kin_R;
            mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin, _n_xi_Kin) =  mat_A1kin * der_kin_R;

            ///// Add the identity matrix to the  the mass and conductance matrix
            mat_p1Ftrans.block(0, 0, _n_xi_Sorp_tilde, _n_xi_Sorp_tilde)                                                      = MathLib::LocalMatrix::Identity(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde);
            mat_p1Ftrans.block(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Min_tilde)                          = MathLib::LocalMatrix::Identity(_n_xi_Min_tilde, _n_xi_Min_tilde);
            mat_p1Ftrans.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde, _n_xi_Sorp)                           = -1.0 * MathLib::LocalMatrix::Identity(_n_xi_Sorp_tilde, _n_xi_Sorp);
            mat_p1Ftrans.block(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min_tilde, _n_xi_Min) = -1.0 * MathLib::LocalMatrix::Identity(_n_xi_Min_tilde, _n_xi_Min);

            //// calculate partial1f. i.e. partial drivative with respect to global variable
            mat_p1F = mat_p1Ftrans - dt * theta_water_content * mat_p1Fder;

            ///// calculate partial2F

            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Sorp, _n_xi_Mob)                          = - dt * theta_water_content * mat_Asorp * der_mob_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, 0, _n_xi_Min, _n_xi_Mob)              = - dt * theta_water_content * mat_Amin * der_mob_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, 0, _n_xi_Kin, _n_xi_Mob)  = - dt * theta_water_content * mat_A1kin * der_mob_R;

            //TODO IMPORTANT: li and ld is merged together. it should be checked later for its correctness.
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Mob, _n_xi_Sorp, _n_xi_Sorp_bar)                         = - dt * theta_water_content * mat_Asorp * der_sorpB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Mob, _n_xi_Min, _n_xi_Sorp_bar)             = - dt * theta_water_content * mat_Amin * der_sorpB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Mob, _n_xi_Kin, _n_xi_Sorp_bar) = - dt * theta_water_content * mat_A1kin * der_sorpB_R;

            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Sorp, _n_xi_Min_bar)                         = - dt * theta_water_content * mat_Asorp * der_minB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min, _n_xi_Min_bar)             = - dt * theta_water_content * mat_Amin * der_minB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Kin, _n_xi_Min_bar) = - dt * theta_water_content * mat_A1kin * der_minB_R;

            // add the relavent identity matrices
            mat_p2F.block(0, _n_xi_Mob, _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_li)                                    = MathLib::LocalMatrix::Identity(_n_xi_Sorp_bar_li, _n_xi_Sorp_bar_li);
            mat_p2F.block(_n_xi_Sorp_bar_li, _n_xi_Mob + _n_xi_Sorp_bar_li, _n_xi_Min_bar, _n_xi_Sorp_bar_ld)    = mat_Ald;
            mat_p2F.block(_n_xi_Sorp_bar_li, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar, _n_xi_Min_bar)           = MathLib::LocalMatrix::Identity(_n_xi_Min_bar, _n_xi_Min_bar);

            ////// calculate vprime
            Vprime(vec_conc, logk_min, mat_S1min, mat_S1mob, mat_S1sorp, mat_S1sorpli, mat_S1kin_ast, mat_S2sorp, mat_vprime);

            // construct local Jacobian matrix
            Jacobian_local = mat_p1F + mat_p2F * mat_vprime;

            // construct global Jacobian matrix
	        std::vector<size_t> node_indx_vec;
	        for(size_t idx = 0; idx != _n_xi_global; idx++)
	        	node_indx_vec[idx]  =  node_idx * idx;

            eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, Jacobian_local);

        }  // end of if statement
    }  // end of node based for loop


    // element based operation: add time and laplas terms
    AddMassLaplasTerms(dt, theta_water_content, eqsJacobian_global);

}

//template <class T1, class T2, class T3>
//void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::NumDiff(std::size_t & col,
//		const double & delta_xi,
//		ogsChem::LocalVector & f,
//		ogsChem::LocalVector & f_old,
//		ogsChem::LocalVector & unknown,
//		ogsChem::LocalVector & DrateDxi)
//		{
//	ogsChem::LocalVector xi = ogsChem::LocalVector::Zero(col);
//	double temp_xi;
//	for(std::size_t i = 0; i < col; i++)
//	{
//		xi        = unknown;
//		temp_xi = xi(i);
//		xi(i)     = xi(i) + (delta_xi * temp_xi);
//		DrateDxi(i) = ( f(xi) - f_old) / (delta_xi * temp_xi);
//	}
//
//		}








template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::Vprime( MathLib::LocalVector & vec_conc,
																  MathLib::LocalVector & logk_min,
																  MathLib::LocalMatrix & mat_S1min,
																  MathLib::LocalMatrix & mat_S1mob,
																  MathLib::LocalMatrix & mat_S1sorp,
																  MathLib::LocalMatrix & mat_S1sorpli,
																  MathLib::LocalMatrix & mat_S1kin_ast,
																  MathLib::LocalMatrix & mat_S2sorp,
																  MathLib::LocalMatrix & mat_vprime)
{
	MathLib::LocalVector  conc_Mob 		  	    = MathLib::LocalVector::Zero(_I_mob);
	MathLib::LocalVector  ln_conc_Mob 		    = MathLib::LocalVector::Zero(_I_mob);
	MathLib::LocalVector  conc_NonMin_bar 		= MathLib::LocalVector::Zero(_I_NMin_bar);
	MathLib::LocalVector  conc_Min_bar	   		= MathLib::LocalVector::Zero(_I_min);
	MathLib::LocalVector  vec_phi				= MathLib::LocalVector::Zero(_n_xi_Min);
	MathLib::LocalMatrix  mat_S1minI, mat_S1minA;
	MathLib::LocalMatrix  mat_A_tilde           = MathLib::LocalMatrix::Zero(_I_mob + _I_NMin_bar, _I_mob + _I_NMin_bar);
	std::size_t i;

	for (i = 0; i <_I_mob + _I_NMin_bar; i++ )
		mat_A_tilde(i,i) = vec_conc(i);

	for (i = 0; i < _I_mob; i++)
	{
		double tmp_x;
		tmp_x    = conc_Mob(i);
		ln_conc_Mob(i)  = std::log(tmp_x);
	}

	vec_phi  = - logk_min + mat_S1min.transpose() * ln_conc_Mob;

	for (i = 0; i < _n_xi_Min; i++)
	{
		//MathLib::LocalVector temp_vec = _mat_S1min.col(i);
		//TODO fix it later
		if(conc_Min_bar(i) >= vec_phi(i)) {
			mat_S1minI.resize( mat_S1min.rows(), mat_S1minI.cols() + 1);
			mat_S1minI.rightCols(1) = mat_S1min.col(i);

		}
		else {
			mat_S1minA.resize( mat_S1min.rows(), mat_S1minA.cols() + 1);
			mat_S1minA.rightCols(1) = mat_S1min.col(i);
		}
	}

	MathLib::LocalMatrix  mat_B  = MathLib::LocalMatrix::Zero(_I_mob + _I_NMin_bar, _n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols());
	mat_B.block(0, 0, _I_mob, _n_xi_Mob)							  =  mat_S1mob;
	mat_B.block(0, _n_xi_Mob, _I_mob, _n_xi_Sorp)					  =  mat_S1sorp;
	mat_B.block(0, _n_xi_Mob + _n_xi_Sorp, _I_mob, mat_S1minI.cols()) =  mat_S1minI;
	mat_B.block(_I_mob, _n_xi_Mob, _I_NMin_bar, _n_xi_Sorp) 		  =  mat_S2sorp;

	MathLib::LocalMatrix  mat_C  = MathLib::LocalMatrix::Zero(_I_mob + _I_NMin_bar, _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin);
	mat_B.block(0, 0, _I_mob, _n_xi_Sorp_tilde)							  =  -1.0 * mat_S1sorpli;
	mat_B.block(0, _n_xi_Sorp_tilde, _I_mob, _n_xi_Min)					  =  -1.0 * mat_S1min;
	mat_B.block(0, _n_xi_Sorp_tilde + _n_xi_Min, _I_mob, _n_xi_Kin)		  =  -1.0 * mat_S1kin_ast;
	MathLib::LocalMatrix J_temp = mat_B.transpose() * mat_A_tilde * mat_B;
	std::size_t sol_size = _n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols();
	//int n = 10000;
	Eigen::SparseMatrix<double> J(sol_size, sol_size);
	//J.reserve(VectorXi::Constant(sol_size,sol_size));

	for(std::size_t i = 0; i < sol_size; i++){
		for(std::size_t j = 0; j < sol_size; j++){
			J.insert(i,j)     = J_temp(i,j);
		}
	}

	Eigen::VectorXd x(sol_size), b(sol_size);
	// solve the linear system
	for(i = 0; i < _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin; i++)
	{
		b = mat_B.transpose() * mat_A_tilde * mat_C.col(i);
		// fill A and b
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
		solver.compute(J);
		x = solver.solve(b);

		mat_vprime.col(i) = x;

	}


}

template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::AddMassLaplasTerms(const double & dt, const double & theta_water_content, LinearSolverType & eqsJacobian_global)
		        {


		    const size_t n_ele = _msh->getNumberOfElements();
		    double _theta(1.0);

		    for (size_t i=0; i<n_ele; i++)
		    {
		        MeshLib::IElement *e = _msh->getElement(i);
		        FemLib::IFiniteElement* _fe;

		        std::vector<size_t> ele_node_ids;
		        //MathLib::Vector<size_t> ele_node_ids;
		        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);

		        //TODO get the water content and multiply it in LHS and RHS

		        //double dt = delta_t.getTimeStepSize();

		        MathLib::LocalMatrix localM = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
		        MathLib::LocalMatrix localK = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
		        MathLib::LocalMatrix localDispersion = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
		        MathLib::LocalMatrix localAdvection = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
		        MathLib::LocalVector F = MathLib::LocalVector::Zero(ele_node_ids.size());
		        MathLib::LocalMatrix poro, d_poro, dispersion_diffusion;
		        MathLib::LocalMatrix v, v2;

		        //assembleODE(time, e, local_u_n1, local_u_n, M, K, F);
		        //_fe = _feObjects->getFeObject(e);

		        const size_t n_dim = e->getDimension();
		        size_t mat_id = e->getGroupID();
		        MaterialLib::PorousMedia* _pm = Ogs6FemData::getInstance()->list_pm[mat_id];

		        localDispersion.setZero(localK.rows(), localK.cols());
		        localAdvection.setZero (localK.rows(), localK.cols());

		        double cmp_mol_diffusion = .0;
		        // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

		        FemLib::IFemNumericalIntegration* _q = _fe->getIntegrationMethod();
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
		            _vel->eval(gp_pos, v);
		            v2 = v.topRows(n_dim).transpose();

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
		            localM(idx_ml, idx_ml) = theta_water_content * mass_lump_val;
		        }



		        std::size_t xi_count, node_indx(0);
		        std::vector<size_t> node_indx_vec;
		        std::vector<size_t>  col_indx_vec;

		        for (xi_count = 0; xi_count < _n_xi_Sorp; xi_count++)
		        {

		        	for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
		        	{
		        		node_indx = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + _n_xi_Min_tilde + xi_count;
		        		std::size_t col_indx  = ele_node_ids[idx] * _n_xi_global + xi_count;
		        		node_indx_vec[idx]    = node_indx;
		        		col_indx_vec[idx]     = col_indx;
		        	}

	        		// add conductance matrix
	        		eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, dt * localK);
	        		// add storage term
	        		eqsJacobian_global.addAsub(node_indx_vec, col_indx_vec, localM);
		        }


		        for (xi_count = 0; xi_count < _n_xi_Min; xi_count++)
		        {

		        	for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
		        	{
		        		node_indx = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + xi_count;
		        		std::size_t col_indx  = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + xi_count;
		        		node_indx_vec[idx]    = node_indx;
		        		col_indx_vec[idx]     = col_indx;
		        	}
	        		// add conductance matrix
	        		eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, dt * localK);
	        		// add storage term
	        		eqsJacobian_global.addAsub(node_indx_vec, col_indx_vec, localM);
		        }


		        for (xi_count = 0; xi_count < _n_xi_Kin; xi_count++)
		        {

		        	for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
		        	{
		        		node_indx = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min + xi_count;
		        		node_indx_vec[idx]    = node_indx;
		        	}

	        		// add mass & conductance matrix
	        		eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, localM + (dt * localK));

		        }

		    }
		}


template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::cal_nodal_rate(ogsChem::LocalVector &xi,
																  ogsChem::LocalVector &local_eta_bar,
																  ogsChem::LocalVector &local_eta,
																  std::size_t & n_xi_Kin_total,
																  ogsChem::LocalMatrix & mat_S1_ast,
																  ogsChem::LocalMatrix & mat_S2_ast,
																  ogsChem::LocalMatrix & mat_S1_orth,
																  ogsChem::LocalMatrix & mat_S2_orth,
																  ogsChem::LocalVector &vec_rate_new )
{
	// declare local temp variable
	ogsChem::LocalVector local_xi, local_xi_bar, local_c_mob, local_c_immob;
	local_xi 	  = ogsChem::LocalVector::Zero(_n_xi_Mob +_n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	local_c_mob   = ogsChem::LocalVector::Zero(_I_mob);
	local_xi_bar  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	local_c_immob = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);


	ogsChem::LocalVector local_xi_global = ogsChem::LocalVector::Zero(_n_xi_global);
	ogsChem::LocalVector local_xi_local  = ogsChem::LocalVector::Zero(_n_xi_local);
	ogsChem::LocalVector local_conc      = ogsChem::LocalVector::Zero(_n_Comp);
    local_xi_global = xi.head(_n_xi_global);
    local_xi_local  = xi.tail(_n_xi_local);

	local_xi.head	(this->_n_xi_Mob) 					  				   = local_xi_local.segment( 0,this->_n_xi_Mob);;
	local_xi.segment(this->_n_xi_Mob, this->_n_xi_Sorp) 		           = local_xi_global.segment(this->_n_xi_Sorp + this->_n_xi_Min,this->_n_xi_Sorp);;
	local_xi.segment(this->_n_xi_Mob + this->_n_xi_Sorp, this->_n_xi_Min ) = local_xi_global.segment(this->_n_xi_Sorp + this->_n_xi_Min + this->_n_xi_Sorp,this->_n_xi_Min);
	local_xi.tail	(this->_n_xi_Kin) 									   = local_xi_global.segment(this->_n_xi_Sorp + this->_n_xi_Min + this->_n_xi_Sorp+ this->_n_xi_Min,this->_n_xi_Kin);

	local_xi_bar.head	(_n_xi_Sorp_bar) 					 =  local_xi_local.segment(this->_n_xi_Mob,this->_n_xi_Sorp_bar);
	local_xi_bar.segment(_n_xi_Sorp_bar, _n_xi_Min_bar) 	 =  local_xi_local.segment(this->_n_xi_Mob+this->_n_xi_Sorp_bar,this->_n_xi_Min_bar);
	local_xi_bar.tail	(_n_xi_Kin_bar) 					 =  local_xi_local.segment(this->_n_xi_Mob+this->_n_xi_Sorp_bar+this->_n_xi_Min_bar,this->_n_xi_Kin_bar);

	local_c_mob   = mat_S1_ast * local_xi   + mat_S1_orth * local_eta;
	local_c_immob = mat_S2_ast * local_xi_bar + mat_S2_orth * local_eta_bar;

	local_conc.topRows( this->_I_mob ) = local_c_mob;
	local_conc.bottomRows( this->_I_NMin_bar + this->_I_min ) = local_c_immob;

//    // testing if the non-negative stablilization will help?
//    for (size_t i=0; i < local_conc.rows(); i++)    {
//        if ( local_conc(i) < 0.0 )
//            local_conc(i) = 1.0e-99;
//    }
//    // end of testing


	// then calculate the rates and fill them in the rate vector
	for (std::size_t i=0; i < n_xi_Kin_total; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( local_conc );
		vec_rate_new(i) = this->_list_kin_reactions[i]->getRate();
	}

}
