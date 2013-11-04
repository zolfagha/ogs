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
#include <fstream>


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
template <class T_DIS_SYS,
          class T_LINEAR_SOLVER,
          class T_USER_FUNCTION_DATA>
class TemplateTransientDxFEMFunction_GIA_Reduct
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;
    typedef typename T_DIS_SYS::template MyLinearEquation<T_LINEAR_SOLVER,DiscreteLib::SparsityBuilderFromNodeConnectivity>::type MyLinaerEQS;
    typedef SolutionLib::FemVariable MyVariable;
    typedef T_LINEAR_SOLVER LinearSolverType;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

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
     	 _n_Comp(_ReductionGIA->get_n_Comp()), _I_mob(_ReductionGIA->get_n_Comp_mob()), _I_min(_ReductionGIA->get_n_Comp_min()), _n_xi_Kin_bar(_ReductionGIA->get_n_xi_Kin_bar()), _I_sorp(_ReductionGIA->get_n_Comp_sorb()),
     	 _n_xi_Mob(_ReductionGIA->get_n_xi_Mob()), _n_eta(_ReductionGIA->get_n_eta()), _n_eta_bar(_ReductionGIA->get_n_eta_bar()), _n_xi_Sorp_tilde(_ReductionGIA->get_n_xi_Sorp_tilde()),
     	_n_xi_Min_tilde(_ReductionGIA->get_n_xi_Min_tilde()), _n_xi_Sorp(_ReductionGIA->get_n_xi_Sorp()), _n_xi_Min(_ReductionGIA->get_n_xi_Min()), _n_xi_Sorp_bar_li(_ReductionGIA->get_n_xi_Sorp_bar_li()),
     	_n_xi_Sorp_bar_ld(_ReductionGIA->get_n_xi_Sorp_bar_ld()), _n_xi_Kin(_ReductionGIA->get_n_xi_Kin()), _n_xi_Min_bar(_ReductionGIA->get_n_xi_Min_bar()),
     	_n_xi_local(_ReductionGIA->get_n_xi_local()), _n_xi_global(_ReductionGIA->get_n_xi_global()),  _J_tot_kin(_ReductionGIA->get_n_xi_Kin_total()), _n_xi_Sorp_bar(_ReductionGIA->get_n_xi_Sorp_bar()),
     	_xi_local_new(userData->get_xi_local_new()), _eta(userData->get_eta()), _eta_bar(userData->get_eta_bar()),
     	_global_vec_Rate(userData->get_global_vec_Rate()), _concentrations(userData->get_concentrations()),  _list_kin_reactions(_ReductionGIA->get_list_kin_reactions())
    {
    };

    ///
    virtual ~TemplateTransientDxFEMFunction_GIA_Reduct() 
    {
        // delete _solv_minimization;
    };

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
    /// @param r       residual
    /// @param du      increment
    void eval(const SolutionLib::SolutionVector &u_n1, const SolutionLib::SolutionVector &r, SolutionLib::SolutionVector &du);

private:
    void GlobalJacobianAssembler(const NumLib::TimeStep & delta_t,
    							 const SolutionLib::SolutionVector & u_cur_xiglob,
    							 LinearSolverType & eqsJacobian_global);
    void Vprime( MathLib::LocalVector & vec_conc,
    			 MathLib::LocalVector & logk_min,
    			 MathLib::LocalMatrix & mat_S1min,
    			 MathLib::LocalMatrix & mat_S1mob,
    			 MathLib::LocalMatrix & mat_S1sorp,
    			 MathLib::LocalMatrix & mat_S1sorpli,
    			 MathLib::LocalMatrix & mat_S1kin_ast,
    			 MathLib::LocalMatrix & mat_S2sorp,
    			 MathLib::LocalMatrix & mat_vprime);
    void NumDiff(std::size_t & col,
    			 const double & delta_xi,
    			 ogsChem::LocalVector & f,
    			 ogsChem::LocalVector & f_old,
    			 ogsChem::LocalVector & unknown,
    			 ogsChem::LocalVector & DrateDxi);
    void AddMassLaplasTerms(const NumLib::TimeStep & delta_t,
    						LinearSolverType & eqsJacobian_global);
	/*
    void cal_nodal_rate(ogsChem::LocalVector &xi,
    					ogsChem::LocalVector &local_eta_bar,
    					ogsChem::LocalVector &local_eta,
    					std::size_t & n_xi_Kin_total,
    					ogsChem::LocalMatrix & mat_S1_ast,
    					ogsChem::LocalMatrix & mat_S2_ast,
    					ogsChem::LocalMatrix & mat_S1_orth,
    					ogsChem::LocalMatrix & mat_S2_orth,
    					ogsChem::LocalVector &vec_rate_new);
	*/
private:
    MeshLib::IMesh* _msh;
    DiscreteLib::DofEquationIdTable* _dofManager;
    MyLinaerEQS* _linear_eqs;
    const NumLib::TimeStep* _t_n1;
    const SolutionLib::SolutionVector* _u_n0;
    std::vector<MyVariable*> _list_var;
    T_USER_FUNCTION_DATA* _userData;
    //MyDiscreteSystem* _dis;
    ogsChem::chemReductionGIA* _ReductionGIA;
    std::vector<MyNodalFunctionScalar*>  &_xi_local_new, &_eta, &_eta_bar, &_global_vec_Rate, &_concentrations;
    std::vector<ogsChem::chemReactionKin*>  _list_kin_reactions;
    FemLib::IFemNumericalIntegration* _q;
    FemLib::IFiniteElement* _fe;
    NumLib::ITXFunction* _vel;
    //SolutionLib::FemDirichletBC *bc1;

    /**
      * pointer to the local problem class.
      */
    LocalProblem* _solv_minimization;

    size_t _n_xi_global, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Sorp, _n_xi_Min, _n_xi_Kin, _n_xi_local, _n_xi_Min_bar, _n_eta, _n_eta_bar, _n_xi_Mob, _n_xi_Kin_bar, _J_tot_kin, _n_xi_Sorp_bar;
    size_t _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_ld, _n_Comp, _I_mob, _I_min, _I_sorp;
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

    LinearSolverType* eqsJacobian_global = _linear_eqs->getLinearSolver();
    // Assembling global jacobian
    GlobalJacobianAssembler(t_n1, u_n1, *eqsJacobian_global);

    // set residual
    _linear_eqs->addRHS(r, -1.0);

    // solve
    _linear_eqs->solve();
    _linear_eqs->getX(du);

}

template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::GlobalJacobianAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, LinearSolverType & eqsJacobian_global)
{

    size_t i, nnodes;
    nnodes = _msh->getNumberOfNodes();
    const double theta_water_content(0.5);
    const double delta_xi = 1E-12;
    std::size_t n_xi_total = _n_xi_local + _n_xi_global;
    // _solv_minimization = new LocalProblem( _ReductionGIA);

    MathLib::LocalVector Unknown_vec 	 = MathLib::LocalVector::Zero(n_xi_total);
    MathLib::LocalMatrix DrateDxi    	 = MathLib::LocalMatrix::Zero(_J_tot_kin, n_xi_total);
    MathLib::LocalVector vec_rate_new    = MathLib::LocalVector::Zero(_J_tot_kin);
    ogsChem::LocalMatrix der_sorpT_R =  MathLib::LocalMatrix::Zero(_J_tot_kin, _n_xi_Sorp_tilde);
    ogsChem::LocalMatrix der_minT_R  =  MathLib::LocalMatrix::Zero(_J_tot_kin, _n_xi_Min_tilde);
    ogsChem::LocalMatrix der_kin_R   =  MathLib::LocalMatrix::Zero(_J_tot_kin, _n_xi_Kin);
    ogsChem::LocalMatrix der_mob_R   = MathLib::LocalMatrix::Zero(_J_tot_kin, _n_xi_Mob);
    ogsChem::LocalMatrix der_sorpB_R = MathLib::LocalMatrix::Zero(_J_tot_kin, _n_xi_Sorp_bar);
    ogsChem::LocalMatrix der_minB_R  = MathLib::LocalMatrix::Zero(_J_tot_kin, _n_xi_Min_bar);

    // current xi global
    MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde,
                         loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;
    // previous xi global
    //MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;

    // current xi local
	MathLib::LocalVector loc_cur_xi_local, loc_cur_xi_Sorp_bar, loc_cur_xi_Min_bar, loc_cur_xi_Sorp_bar_li, loc_cur_xi_Sorp_bar_ld, tmp_xi_global, tmp_xi_local, tmp_vec_conc;

    // eta vectors
    MathLib::LocalVector loc_cur_eta, loc_cur_eta_bar;

    // local concentration vector
    MathLib::LocalVector vec_conc;

    MathLib::LocalMatrix mat_p1Fder,mat_p1Ftrans, mat_p1F, mat_p2F, mat_Global_Jacobian;
    //std::size_t n_xi_Kin_total = _ReductionGIA->get_n_xi_Kin_total();

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
    // current eta mobile and immobile
    loc_cur_eta             = MathLib::LocalVector::Zero( _n_eta );
    loc_cur_eta_bar         = MathLib::LocalVector::Zero( _n_eta_bar );
    vec_conc                = MathLib::LocalVector::Zero(_n_Comp);

    MathLib::LocalMatrix mat_vprime 	 = MathLib::LocalMatrix::Zero(_n_xi_local,_n_xi_global);
    MathLib::LocalMatrix Jacobian_local  = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    MathLib::LocalVector vec_rate_old 	 = MathLib::LocalVector::Zero(_J_tot_kin);
    MathLib::LocalVector vec_Rate 		 = MathLib::LocalVector::Zero(_J_tot_kin);
    MathLib::LocalVector logk_min	   	 = MathLib::LocalVector::Zero(_n_xi_Min);
    MathLib::LocalMatrix mat_S1min, mat_S1mob, mat_S1sorp, mat_S1sorpli, mat_S1kin_ast, mat_S2sorp;

    mat_p1Fder  = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p1Ftrans = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p1F     = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p2F     = MathLib::LocalMatrix::Zero(_n_xi_global, _n_xi_local);
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

    std::vector<size_t> node_indx_vec;
    for (i=0; i< _n_xi_global; i++)
        node_indx_vec.push_back(0); 

    // loop over all the nodes
    for (size_t node_idx=0; node_idx < nnodes; node_idx++ )
    {
            // on each node, get the right start value
            // get the right set of xis
            for (i=0; i < _n_xi_global; i++)
                loc_cur_xi_global[i] = u_cur_xiglob[node_idx * _n_xi_global + i];
            for (i=0; i < _n_xi_local; i++)
                loc_cur_xi_local[i] = _xi_local_new[i]->getValue(node_idx);
            for (i=0; i < _n_eta; i++)
                loc_cur_eta[i] = _eta[i]->getValue(node_idx);
            for (i=0; i < _n_eta_bar; i++)
                loc_cur_eta_bar[i] = _eta_bar[i]->getValue(node_idx);

            for (i=0; i < _J_tot_kin; i++)
                vec_rate_old[i] = _global_vec_Rate[i]->getValue(node_idx);

            for (i=0; i < _n_Comp; i++)
                vec_conc[i] = this->_concentrations[i]->getValue(node_idx);

			tmp_vec_conc = vec_conc; // just initialize memory. 

            double DrateDxi_temp (0.0);
            Unknown_vec.head(_n_xi_global) = loc_cur_xi_global;
            Unknown_vec.tail(_n_xi_local) = loc_cur_xi_local;

            if(_n_xi_Kin != 0)
            {
            	for(std::size_t i = 0; i < n_xi_total ; i++ )
            	{
            		ogsChem::LocalVector xi = ogsChem::LocalVector::Zero(n_xi_total);
            		xi        = Unknown_vec;

            		if (Unknown_vec.norm() >= 1.0e-16)
            			xi(i)     += delta_xi * Unknown_vec.norm();
            		else
            			xi(i)     += delta_xi;

					tmp_xi_global = xi.head(_n_xi_global); 
					tmp_xi_local = xi.tail(_n_xi_local);
					// convert xi and eta values to vec_conc
					_ReductionGIA->EtaXi2Conc(loc_cur_eta, loc_cur_eta_bar, tmp_xi_global, tmp_xi_local, tmp_vec_conc); 
            		// calculate new rate term
					_ReductionGIA->Calc_Kin_Rate(tmp_vec_conc,
												 vec_rate_new);

            		for(std::size_t j = 0; j < _J_tot_kin ; j++ )
            		{
            			if (std::abs(Unknown_vec(i)) >= 1.0e-16)
            				DrateDxi_temp = ( vec_rate_new (j)- vec_rate_old(j)) / delta_xi / Unknown_vec.norm();
            			else
            				DrateDxi_temp = ( vec_rate_new (j)- vec_rate_old(j)) / delta_xi;

            			DrateDxi(j,i)  = DrateDxi_temp;
            		}
            	}

            	der_sorpT_R 					 =  DrateDxi.block(0, 0, _J_tot_kin, _n_xi_Sorp_tilde);
            	der_minT_R 				    	 =  DrateDxi.block(0, _n_xi_Sorp_tilde, _J_tot_kin, _n_xi_Min_tilde);
            	der_kin_R 				    	 =  DrateDxi.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _J_tot_kin, _n_xi_Kin);  // try der_kin_R in the laplas assembly

            	der_mob_R 				    	 =  DrateDxi.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin,_J_tot_kin,_n_xi_Mob);
            	der_sorpB_R 				     =  DrateDxi.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin + _n_xi_Mob, _J_tot_kin, _n_xi_Sorp_bar);
            	der_minB_R	 				     =  DrateDxi.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin + _n_xi_Mob + _n_xi_Sorp_bar, _J_tot_kin, _n_xi_Min_bar);



            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Sorp, _n_xi_Sorp_tilde)                         =  mat_Asorp * der_sorpT_R;
            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, 0, _n_xi_Min, _n_xi_Sorp_tilde)             =  mat_Amin * der_sorpT_R;
            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, 0, _n_xi_Kin, _n_xi_Sorp_tilde) =  mat_A1kin * der_sorpT_R;

            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde, _n_xi_Sorp, _n_xi_Min_tilde)                         =  mat_Asorp * der_minT_R;
            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Sorp_tilde, _n_xi_Min, _n_xi_Min_tilde)             =  mat_Amin * der_minT_R;
            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp_tilde, _n_xi_Kin, _n_xi_Min_tilde) =  mat_A1kin * der_minT_R;

            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp, _n_xi_Kin)                         =  mat_Asorp * der_kin_R;
            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Min, _n_xi_Kin)             =  mat_Amin * der_kin_R;
            	mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin, _n_xi_Kin)   =  mat_A1kin * der_kin_R;

			} // end of if (_n_xi_Kin != 0)

            // Add the identity matrix to the  the mass and conductance matrix
            mat_p1Ftrans.block(0, 0, _n_xi_Sorp_tilde, _n_xi_Sorp_tilde)                                                      = MathLib::LocalMatrix::Identity(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde);
            mat_p1Ftrans.block(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Min_tilde)                          = MathLib::LocalMatrix::Identity(_n_xi_Min_tilde, _n_xi_Min_tilde);
            mat_p1Ftrans.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde, _n_xi_Sorp)                           = -1.0 * MathLib::LocalMatrix::Identity(_n_xi_Sorp_tilde, _n_xi_Sorp);
            mat_p1Ftrans.block(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min_tilde, _n_xi_Min) = -1.0 * MathLib::LocalMatrix::Identity(_n_xi_Min_tilde, _n_xi_Min);

            // debugging--------------------------
            // std::cout << "======================================== \n";
            // std::cout << "mat_p1Ftrans: \n";
            // std::cout << mat_p1Ftrans << std::endl;
            // std::cout << "======================================== \n";
            // end of debugging-------------------

            //// calculate partial1f. i.e. partial drivative with respect to global variable
			mat_p1F = mat_p1Ftrans  - theta_water_content * mat_p1Fder;  // dt is mulitpiled by mass matrix

            // debugging--------------------------
            // std::cout << "======================================== \n";
            // std::cout << "mat_p1F: \n";
            // std::cout << mat_p1F << std::endl;
            // std::cout << "======================================== \n";
            // end of debugging-------------------

            // calculate partial2F
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Sorp, _n_xi_Mob)                          = - theta_water_content * mat_Asorp * der_mob_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, 0, _n_xi_Min, _n_xi_Mob)              = - theta_water_content * mat_Amin * der_mob_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, 0, _n_xi_Kin, _n_xi_Mob)  = - theta_water_content * mat_A1kin * der_mob_R;

            //TODO IMPORTANT: li and ld is merged together. it should be checked later for its correctness.
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Mob, _n_xi_Sorp, _n_xi_Sorp_bar)                         = - theta_water_content * mat_Asorp * der_sorpB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Mob, _n_xi_Min, _n_xi_Sorp_bar)             = - theta_water_content * mat_Amin * der_sorpB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Mob, _n_xi_Kin, _n_xi_Sorp_bar) = - theta_water_content * mat_A1kin * der_sorpB_R;

            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Sorp, _n_xi_Min_bar)                         = - theta_water_content * mat_Asorp * der_minB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min, _n_xi_Min_bar)             = - theta_water_content * mat_Amin * der_minB_R;
            mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Kin, _n_xi_Min_bar) = - theta_water_content * mat_A1kin * der_minB_R;

            // debugging--------------------------
            // std::cout << "======================================== \n";
            // std::cout << "mat_p2F: \n";
            // std::cout << mat_p2F << std::endl;
            // std::cout << "======================================== \n";
            // end of debugging-------------------

            // add the relavent identity matrices
            mat_p2F.block(0, _n_xi_Mob, _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_li)                                    = MathLib::LocalMatrix::Identity(_n_xi_Sorp_bar_li, _n_xi_Sorp_bar_li);
            mat_p2F.block(_n_xi_Sorp_bar_li, _n_xi_Mob + _n_xi_Sorp_bar_li, _n_xi_Min_bar, _n_xi_Sorp_bar_ld)    = mat_Ald;
            mat_p2F.block(_n_xi_Sorp_bar_li, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar, _n_xi_Min_bar)           = MathLib::LocalMatrix::Identity(_n_xi_Min_bar, _n_xi_Min_bar);

            // calculate vprime
            Vprime(vec_conc, logk_min, mat_S1min, mat_S1mob, mat_S1sorp, mat_S1sorpli, mat_S1kin_ast, mat_S2sorp, mat_vprime);

            // construct local Jacobian matrix
            Jacobian_local = mat_p1F + mat_p2F * mat_vprime;

            // debugging--------------------------
            // std::cout << "======================================== \n";
            // std::cout << "mat_vprime: \n";
            // std::cout << mat_vprime << std::endl;
            // std::cout << "Jacobian_local: \n";
            // std::cout << Jacobian_local << std::endl;
            // std::cout << "======================================== \n";
            // end of debugging-------------------

            // construct global Jacobian matrix
            for(size_t idx = 0; idx < _n_xi_global; idx++)
            {
            	node_indx_vec[idx] = node_idx * _n_xi_global + idx;
            }

            eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, Jacobian_local);
            mat_vprime.setZero();
			mat_p2F.setZero(); 
			mat_p1F.setZero();
			mat_p1Ftrans.setZero(); 
			mat_p1Fder.setZero(); 

            }  // end of node based for loop
            node_indx_vec.clear();

            //// --------debugging--------------
            // std::ofstream globalJ_node ("globalJ_node.txt");
            // eqsJacobian_global.printout(globalJ_node);
            //// --------end of debugging-------

    // element based operation: add time and laplas terms
    AddMassLaplasTerms(delta_t, eqsJacobian_global);

    //// --------debugging--------------
    // std::ofstream globalJ ("globalJ.txt");
    // eqsJacobian_global.printout(globalJ);
    //// --------end of debugging-------


}


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
	MathLib::LocalVector  ln_conc_Mob 		    = MathLib::LocalVector::Zero(_I_mob);
	MathLib::LocalVector  conc_Min_bar	   		= MathLib::LocalVector::Zero(_I_min);
	MathLib::LocalVector  vec_phi				= MathLib::LocalVector::Zero(_n_xi_Min);
	MathLib::LocalMatrix  mat_S1minI            = MathLib::LocalMatrix::Zero(_I_mob,0);
	MathLib::LocalMatrix  mat_S1minA            = MathLib::LocalMatrix::Zero(_I_mob,0);
	MathLib::LocalMatrix  mat_A_tilde           = MathLib::LocalMatrix::Zero(_I_mob + _I_sorp, _I_mob + _I_sorp);
	std::size_t i;

	for (i = 0; i < _I_mob + _I_sorp; i++)
		mat_A_tilde(i,i) = 1.0 / vec_conc(i);

	for (i = 0; i < _I_mob; i++)
	{
		double tmp_x;
		tmp_x    = vec_conc(i);
		ln_conc_Mob(i)  = std::log(tmp_x);
	}

	vec_phi		  = - logk_min + mat_S1min.transpose() * ln_conc_Mob;
	conc_Min_bar  = vec_conc.tail(_I_min);

    mat_S1minI.setZero();
    mat_S1minA.setZero();

	for (i = 0; i < _n_xi_Min; i++)
	{
		if(conc_Min_bar(i) >= vec_phi(i)) {
			mat_S1minI.resize( mat_S1min.rows(), mat_S1minI.cols() + 1);
			mat_S1minI.rightCols(1) = mat_S1min.col(i);
		}
		else {
			mat_S1minA.resize( mat_S1min.rows(), mat_S1minA.cols() + 1);
			mat_S1minA.rightCols(1) = mat_S1min.col(i);
		}
	}

	MathLib::LocalMatrix  mat_B = MathLib::LocalMatrix::Zero(_I_mob + _I_sorp, _n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols());
	mat_B.block(0, 0, _I_mob, _n_xi_Mob)							  =  mat_S1mob;
	mat_B.block(0, _n_xi_Mob, _I_mob, _n_xi_Sorp)					  =  mat_S1sorp;
	mat_B.block(0, _n_xi_Mob + _n_xi_Sorp, mat_S1minI.rows(), mat_S1minI.cols()) =  mat_S1minI;
	mat_B.block(_I_mob, _n_xi_Mob, _I_sorp, _n_xi_Sorp) 		  =  mat_S2sorp;

	MathLib::LocalMatrix  mat_C = MathLib::LocalMatrix::Zero(_I_mob + _I_sorp, _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin);
	mat_C.block(0, 0, _I_mob, _n_xi_Sorp_tilde)							  =  - 1.0 * mat_S1sorpli;
	mat_C.block(0, _n_xi_Sorp_tilde, _I_mob, _n_xi_Min)					  =  - 1.0 * mat_S1min;
	mat_C.block(0, _n_xi_Sorp_tilde + _n_xi_Min, _I_mob, _n_xi_Kin)		  =  - 1.0 * mat_S1kin_ast;

	MathLib::LocalMatrix J_temp = MathLib::LocalMatrix::Zero( _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin,  _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin);
	J_temp					    = mat_B.transpose() * mat_A_tilde * mat_B;

	std::size_t sol_size 		= _n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols();

    // using minimization solver of the local problem.

    MathLib::LocalVector b 	= MathLib::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols());
    //MathLib::LocalVector dx = MathLib::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols());
    MathLib::LocalMatrix x  = MathLib::LocalMatrix::Zero(_n_xi_Mob + _n_xi_Sorp + mat_S1minI.cols(), _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin);
    std::size_t indx_dummy = 0;
    // solve the linear system
    for(i = 0; i < _n_xi_Sorp_tilde + _n_xi_Min + _n_xi_Kin; i++)
    {
	    b = mat_B.transpose() * mat_A_tilde * mat_C.col(i);

	    // _solv_minimization->Solv_Minimization(indx_dummy, J_temp, b, dx);  // b is multiplied by -1 in this function.
	    // x.col(i) = dx;

	    // using the standard direct solver
	    x.col(i) = J_temp.fullPivHouseholderQr().solve( b );

    }

	MathLib::LocalMatrix cols_xi_sorp_tilde   = x.block(0, 0, sol_size, _n_xi_Sorp_tilde);
	mat_vprime.block(0, 0, sol_size, _n_xi_Sorp_tilde) = cols_xi_sorp_tilde;

	MathLib::LocalMatrix cols_xi_min_tilde   = x.block(0, _n_xi_Sorp_tilde, sol_size, _n_xi_Min_tilde);
	mat_vprime.block(0, _n_xi_Sorp_tilde, sol_size, _n_xi_Min_tilde ) = cols_xi_min_tilde;

	//MathLib::LocalMatrix cols_xi_kin  = x.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde, sol_size, _n_xi_Kin);  //RZ: wrong, matrix x contains xi mob,sorp,minI
	MathLib::LocalMatrix cols_xi_kin  = MathLib::LocalMatrix::Identity(_n_xi_Kin_bar, _n_xi_Kin_bar);
	//mat_vprime.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, sol_size, _n_xi_Kin) = cols_xi_kin;
	mat_vprime.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin_bar, _n_xi_Kin_bar) = cols_xi_kin;

	mat_B.setZero();
}

template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>
     ::AddMassLaplasTerms(const NumLib::TimeStep & delta_t, 
                          LinearSolverType & eqsJacobian_global)
{
    const size_t n_ele = _msh->getNumberOfElements();
    std::size_t n_dim, mat_id; 
    std::size_t i, j, xi_count, idx;
    int idx_ml; 
    double _theta(1.0);
    double dt = delta_t.getTimeStepSize();
    double cmp_mol_diffusion;
    MeshLib::IElement *e; 
    MaterialLib::PorousMedia* _pm;
    std::vector<size_t> ele_node_ids;
    std::vector<size_t> node_indx_vec;
	std::vector<size_t>  col_indx_vec;

    MathLib::LocalMatrix localM, localK, localDispersion, localAdvection; 
    MathLib::LocalVector F; 
    MathLib::LocalMatrix dispersion_diffusion, d_poro, poro; 
    NumLib::ITXFunction::DataType v;

    for ( i=0; i<n_ele; i++)
    {
        e = _msh->getElement(i);
		e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);

	    localM          = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
	    localK          = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
	    localDispersion = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
	    localAdvection  = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
	    F               = MathLib::LocalVector::Zero(ele_node_ids.size());
	    d_poro          = MathLib::LocalMatrix::Zero(3,3);
        poro            = MathLib::LocalMatrix::Zero(1,1);

	    //assembleODE(time, e, local_u_n1, local_u_n, M, K, F);
	    _fe = _userData->get_feObjects()->getFeObject(*e);

	    n_dim = e->getDimension();
	    mat_id = e->getGroupID();
	    _pm = Ogs6FemData::getInstance()->list_pm[mat_id];

	    localDispersion.setZero(localK.rows(), localK.cols());
	    localAdvection.setZero (localK.rows(), localK.cols());

	    cmp_mol_diffusion = .0;
	    // _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

	    _q = _fe->getIntegrationMethod();
	    double gp_x[3], real_x[3];
	    poro.setZero();
	    d_poro.setZero();
	    double disp_l = 0.0;
	    double disp_t = 0.0;

		for ( j=0; j < _q->getNumberOfSamplingPoints(); j++)
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
		    dispersion_diffusion *= disp_t * v.norm();
		    dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm();
		    dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);

		    _fe->integrateWxN(j, poro, localM);
		    _fe->integrateDWxDN(j, dispersion_diffusion, localDispersion);
		    _fe->integrateWxDN(j, v2, localAdvection);
		}  // end of for loop over j

		localK = localDispersion + localAdvection;
		
		// disable mass lumping
		/*
		// mass lumping----------------------------
		for ( idx_ml=0; idx_ml < localM.rows(); idx_ml++ )
		{
		    double mass_lump_val;
		    mass_lump_val = localM.row(idx_ml).sum();
		    localM.row(idx_ml).setZero();
		    localM(idx_ml, idx_ml) =  mass_lump_val;
		}
		*/

		node_indx_vec.resize(ele_node_ids.size());
		col_indx_vec.resize(ele_node_ids.size());
                
		for (xi_count = 0; xi_count < _n_xi_Sorp; xi_count++)
		{

		    for( idx = 0; idx < ele_node_ids.size(); idx++)
		    {
		        node_indx_vec[idx] = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + _n_xi_Min_tilde + xi_count;
		        col_indx_vec[idx]   = ele_node_ids[idx] * _n_xi_global + xi_count;
		    }
	        // add conductance matrix
	        eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, localK );
	        // add storage term
	        eqsJacobian_global.addAsub(node_indx_vec, col_indx_vec, 1.0 / dt * localM);
		}


		for (xi_count = 0; xi_count < _n_xi_Min; xi_count++)
		{

		    for( idx = 0; idx < ele_node_ids.size(); idx++)
		    {
		        node_indx_vec[idx] = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + xi_count;
		        col_indx_vec[idx]  = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + xi_count;
		    }
	        eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, localK );
	        // add storage term
	        eqsJacobian_global.addAsub(node_indx_vec, col_indx_vec, 1.0 / dt * localM);
		}

		for (xi_count = 0; xi_count < _n_xi_Kin; xi_count++)
		{

		    for( idx = 0; idx < ele_node_ids.size(); idx++)
		    {
		        node_indx_vec[idx]  = ele_node_ids[idx] * _n_xi_global + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min + xi_count;
		    }
	        // add mass & conductance matrix
	        eqsJacobian_global.addAsub(node_indx_vec, node_indx_vec, 1.0 / dt * localM + localK );

		}
		node_indx_vec.clear();
		col_indx_vec.clear();

	}  // end of for i over all elements
}  // end of function AddMassLaplasTerms

/*
template <class T1, class T2, class T3>
void TemplateTransientDxFEMFunction_GIA_Reduct<T1,T2,T3>::cal_nodal_rate(ogsChem::LocalVector &xi,
																  ogsChem::LocalVector &local_eta_bar,
																  ogsChem::LocalVector &local_eta,
																  ogsChem::LocalVector &vec_rate_new )
{
    // declare local temp variable
	ogsChem::LocalVector local_xi_global = xi.head(_n_xi_global);
	ogsChem::LocalVector local_xi_local = xi.tail(_n_xi_local);

	_ReductionGIA->Calc_Kin_Rate(local_xi_local,
		                         local_xi_global,
								 local_eta,
								 local_eta_bar,
								 vec_rate_new);

}
*/

