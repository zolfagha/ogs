/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemEqReactSys.cpp
 *
 * Created on 2013-05-26 by Reza Zolfaghari & Haibing Shao
 */

#include "LocalProblem.h"
#include "logog.hpp"
#define _DEBUG

LocalProblem::LocalProblem(ogsChem::chemReductionGIA* ReductionGIA)
	: _n_Comp(ReductionGIA->get_n_Comp()), _ReductionGIA(ReductionGIA), _I_mob(_ReductionGIA->get_n_Comp_mob()), _I_sorp(_ReductionGIA->get_n_Comp_sorb()), _I_min(_ReductionGIA->get_n_Comp_min()), _n_xi_Kin_bar(_ReductionGIA->get_n_xi_Kin_bar())
	, _n_xi_Mob(_ReductionGIA->get_n_xi_Mob()), _n_eta(_ReductionGIA->get_n_eta()), _n_eta_bar(_ReductionGIA->get_n_eta_bar()), _n_xi_Sorp_tilde(_ReductionGIA->get_n_xi_Sorp_tilde()), _n_xi_Min_tilde(_ReductionGIA->get_n_xi_Min_tilde())
    , _n_xi_Sorp(_ReductionGIA->get_n_xi_Sorp()), _n_xi_Min(_ReductionGIA->get_n_xi_Min()), _n_xi_Sorp_bar_li(_ReductionGIA->get_n_xi_Sorp_bar_li()), _n_xi_Sorp_bar_ld(_ReductionGIA->get_n_xi_Sorp_bar_ld()), _n_xi_Kin(_ReductionGIA->get_n_xi_Kin())
	, _mat_c_mob_2_xi_mob(_ReductionGIA->get_matrix_C2Xi()), _mat_c_immob_2_xi_immob(_ReductionGIA->get_matrix_Cbar2XiBar()), _n_xi_Sorp_bar(_ReductionGIA->get_n_xi_Sorp_bar()), _mat_S1min (_ReductionGIA->get_matrix_S1min())
    , _mat_Ald(_ReductionGIA->get_matrix_Ald()),_mat_Ssorp(_ReductionGIA->get_matrix_Ssorp()), _mat_A2kin(_ReductionGIA->get_matrix_A2kin()), _mat_S1mob(_ReductionGIA->get_matrix_S1mob()),_mat_c_mob_2_eta_mob(_ReductionGIA->get_matrix_C2Eta())
    , _logk_mob(_ReductionGIA->get_logk_mob()), _logk_sorp(_ReductionGIA->get_logk_sorp()), _logk_min(_ReductionGIA->get_logk_min()), _n_xi_Min_bar(_ReductionGIA->get_n_xi_Min_bar()), _I_NMin_bar(_ReductionGIA->get_n_Comp_NMin_bar())
    , _mat_c_immob_2_eta_immob(_ReductionGIA->get_matrix_C2EtaBar())

{

}

LocalProblem::~LocalProblem()
{}


void LocalProblem::solve_LocalProblem_Newton_LineSearch(ogsChem::LocalVector & vec_conc,
														ogsChem::LocalVector & vec_tot_mass_constrain,
														std::size_t & node_idx ,
														double deltaT,
														const double iter_tol,
														const double rel_tol,
														const double max_iter)
{
    ogsChem::LocalVector x, x_new, vec_residual, vec_AI;
    ogsChem::LocalVector dx;
    ogsChem::LocalVector conc_Mob, conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar;
//TODO initialize water content and deltaT
    // number of iterations
    size_t j, iter, n_unknowns;
    const double alpha (0.5);
    double d_norm, d1_norm;

    n_unknowns     = _n_Comp + _n_xi_Kin_bar;
    x              = ogsChem::LocalVector::Ones( n_unknowns );
    x_new          = ogsChem::LocalVector::Ones( n_unknowns );
    dx             = ogsChem::LocalVector::Ones( n_unknowns );
    _mat_Jacobian  = ogsChem::LocalMatrix::Zero(n_unknowns, n_unknowns);
    vec_residual   = ogsChem::LocalVector::Ones( n_unknowns );
    // initialize the _AI vector
    vec_AI		   = ogsChem::LocalVector::Zero( _I_min );

	// vec_tot_mass_constrain contains xi global and eta mobile which acts as a total mass constrain for the local problem.
	_vec_eta			= vec_tot_mass_constrain.segment(0,_n_eta);
	_vec_XiSorpTilde    = vec_tot_mass_constrain.segment(_n_eta,_n_xi_Sorp_tilde);
	_vec_XiMinTilde     = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde,_n_xi_Min_tilde);
	_vec_Xikin          = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Kin);
	_vec_etabar 		= vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin,_n_eta_bar);
	_vec_XiBarKin_old   = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_eta_bar, _n_xi_Kin_bar);

    //concentration vector components, ie, mobile, sorbed, mineral concentrations
    conc_Mob     	 = vec_conc.head(_I_mob );
    //NonMin refers to sorbed concentrations.
    conc_NonMin_bar  = vec_conc.segment(_I_mob, _I_NMin_bar );
    conc_Min_bar     = vec_conc.segment(_I_mob + _I_NMin_bar, _I_min );
    // xi kin bar is both unknown and mass constrain!
    Xi_Kin_bar       = vec_conc.tail(_n_xi_Kin_bar);


    #ifdef _DEBUG
        /*
        std::cout << "Total mass for the basis: \n";
        std::cout << total_mass << "\n"; 
        std::cout << "Conc_basis: \n";
        std::cout << conc_basis << "\n"; 
        std::cout << "Conc_second: \n";
        std::cout << conc_second << std::endl; 
        */
    #endif



    // start solving the system
    iter = 0; 

    // now updating the saturation index and minerals
    this->update_minerals_conc_AI( x, vec_AI );

    // evaluate the residual
    this->calc_residual(x, vec_AI, vec_residual);

    // evaluate norm of residual vector
    d_norm = vec_residual.norm();

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Residual Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif


    while (true)
    {
        #ifdef _DEBUG
            // display the residual
             std::cout << "Iteration #" << iter << "||res|| = " << vec_residual.norm() << "||delta_x|| = " << dx.norm() << std::endl;
        #endif
        // convergence criteria
        if ( d_norm < iter_tol )
        {
            #ifdef _DEBUG
                 std::cout << "Newton iteration successfully converged!\n";
            #endif

            break;  // break the loop
        }
        else if ( dx.norm() < rel_tol )
        {
            #ifdef _DEBUG
            std::cout << "Warning, Newton iteration stagnent on Node #" << node_idx << "! Exit the iteration!\n" ; 
            #endif

            break;  // break the loop
        }
        else if ( iter > max_iter )
        {
            #ifdef _DEBUG
            std::cout << "ERROR! Node #" << node_idx  << "Newton iterationan does not converge! Simulation stops!\n"; 
            #endif

            return; // stop the program
        }
        // form Jacobian matrix
        this->calc_Jacobian(x, vec_AI, vec_residual);
        // solving for increment
        this->Solv_Minimization( node_idx, _mat_Jacobian, vec_residual, dx );

        // increment of unkowns
        this->increment_unknown( x, dx, x_new ); 

        // evaluate residual with x_new
        this->calc_residual(x_new, vec_AI, vec_residual);

        // line search begins
        j = 0; 
        while ( j < max_iter )
        {
            // d1_norm = norm(res,inf);
            d1_norm = vec_residual.norm();
            
            if (d1_norm < d_norm)
                break;
            
            // updating dx
            dx = dx * alpha;
            // increment of unknowns
            this->increment_unknown( x, dx, x_new ); 
            // now updating the saturation index and minerals
            this->update_minerals_conc_AI( x_new, vec_AI );
            // evaluate residual with x_new
            this->calc_residual(x_new, vec_AI, vec_residual);

		#ifdef _DEBUG
            // display the residual
            std::cout << "Line Search Iteration #" << iter << "||res|| = " << vec_residual.norm() << "||delta_x|| = " << dx.norm() << std::endl;
		#endif

            j++;
        }  // end of while
        d_norm = d1_norm; 
        x = x_new; 

        // increase the iteration count
        iter++; 
    }
}

void LocalProblem::calc_residual(ogsChem::LocalVector & vec_unknowns,
								 ogsChem::LocalVector & vec_AI,
								 ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector conc_Mob, conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar;

	//


//
//
    conc_Mob         = vec_unknowns.segment(0, _I_mob );
    //NonMin refers to sorbed concentrations.
    conc_NonMin_bar  = vec_unknowns.segment(_I_mob, _I_NMin_bar );
    conc_Min_bar     = vec_unknowns.segment(_I_mob + _I_NMin_bar, _I_min );
    Xi_Kin_bar   	 = vec_unknowns.segment(_I_mob + _I_NMin_bar + _I_min, _n_xi_Kin_bar );

    // Eq. 3.55
    this->residual_conc_Mob			(conc_Mob, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_conc_Mob Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.56
    this->residual_Eta				(conc_Mob, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_Eta Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.57 - 58
    this->residual_xi_Sorp_tilde	(conc_Mob, conc_NonMin_bar, conc_Min_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_xi_Sorp_tilde Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.59
    this->residual_xi_Min_tilde		(conc_Mob, conc_NonMin_bar, conc_Min_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_xi_Min_tilde Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.60
    this->residual_xi_Kin			(conc_Mob, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_xi_Kin Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.61
    this->residual_conc_Sorp		(conc_Mob, conc_NonMin_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_conc_Sorp Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.62
    this->residual_conc_Min			(conc_Mob, vec_AI, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_conc_Min Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.63
    this->residual_Eta_bar			(conc_NonMin_bar, conc_Min_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_Eta_bar Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.64
    this->residual_xi_KinBar_Eq		(conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_xi_KinBar_Eq Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.65
    this->residual_xi_KinBar_Kin	(conc_Mob, conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "residual_xi_KinBar_Kin Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

}  // end of function calc_residual

void LocalProblem::calc_Jacobian(ogsChem::LocalVector & vec_x,
							   ogsChem::LocalVector & vec_AI,
							   ogsChem::LocalVector & vec_residual)
{
	const double delta_xi = 1.0e-6;
    std::size_t i;
    ogsChem::LocalVector vec_x_incremented, vec_residual_incremented;
    vec_residual_incremented = vec_residual;

	for (i = 0; i < _mat_Jacobian.cols(); i++)
	{
		vec_x_incremented  = vec_x;

//		// numerical protection
//		if( vec_x_incremented.norm() < 1.0e-16)
//		{
//			vec_x_incremented(i) += delta_xi * vec_x_incremented.norm();
//			this->calc_residual(vec_x_incremented,vec_AI,vec_residual_incremented);
//			_mat_Jacobian.col(i) = (vec_residual_incremented - vec_residual ) / (delta_xi * vec_x_incremented.norm());
//
//		}
//		else
		{
		vec_x_incremented(i) += delta_xi;

		this->calc_residual(vec_x_incremented,vec_AI,vec_residual_incremented);
		_mat_Jacobian.col(i) = (vec_residual_incremented - vec_residual ) / delta_xi;
		}
	}

//	using namespace std::placeholders; //for _1, _2, _3...

//	_vec_eta			= vec_tot_mass_constrain.segment(0,_n_eta);
//	_vec_XiSorpTilde    = vec_tot_mass_constrain.segment(_n_eta,_n_xi_Sorp_tilde);
//	_vec_XiMinTilde     = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde,_n_xi_Min_tilde);
//	_vec_Xikin          = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Kin);
//	_vec_etabar 		= vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin,_n_eta_bar);
//	_vec_XiBarKin  	    = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_eta_bar, _n_xi_Kin_bar);
//	_vec_XiBarKin_old   = vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_eta_bar + _n_xi_Kin_bar, _n_xi_Kin_bar);
//
//    _conc_Mob     	 = vec_unknowns.segment(0, _I_mob );
//    //NonMin refers to sorbed concentrations.
//    _conc_NonMin_bar  = vec_unknowns.segment(_I_mob, _I_NMin_bar );
//    _conc_Min_bar     = vec_unknowns.segment(_I_mob + _I_NMin_bar, _I_min );
//    _Xi_Kin_bar   	 = vec_unknowns.segment(_I_mob + _I_NMin_bar + _I_min, _n_xi_Kin_bar );





    //TODO put them on the right rows of the jacobian matrix

//    // Eq. 3.55
//    _mat_Jacobi.block(0, 0, _n_xi_Mob,_I_mob)
//    this->Num_Diff(_n_xi_Mob, _I_mob, _conc_Mob, delta_xi, std::bind(residual_conc_Mob));
//
//    //Eq. 3.56
//    _mat_Jacobi.block(_n_xi_Mob, 0, _n_eta, _I_mob)
//    this->Num_Diff(_n_eta,_I_mob, _conc_Mob, delta_xi,std::bind(residual_conc_Mob,_2, _vec_eta));
//
//    //Eq. 3.57-8
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta, 0, _n_xi_Sorp_tilde, _I_mob)
//	this->Num_Diff(_n_xi_Sorp,_I_mob,_conc_Mob, delta_xi,std::bind(residual_xi_Sorp_tilde,_2, _3, _4, _conc_NonMin_bar, conc_Min_bar, _vec_XiSorpTilde));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta, _I_mob, _n_xi_Sorp_tilde, _I_NMin_bar)
//	this->Num_Diff(_n_xi_Sorp,_I_NMin_bar,_conc_NonMin_bar, delta_xi,std::bind(residual_xi_Sorp_tilde,_1, _3, _4, _conc_Mob, conc_Min_bar, _vec_XiSorpTilde));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta, _I_mob + _I_NMin_bar, _n_xi_Sorp_tilde, _I_min)
//	this->Num_Diff(_n_xi_Sorp,_I_min,_conc_Min_bar, delta_xi,std::bind(residual_xi_Sorp_tilde,_1, _2, _4, _conc_Mob, _conc_NonMin_bar, _vec_XiSorpTilde));
//
//	// Eq.3.59
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde, 0, _n_xi_Min_tilde, _I_mob)
//	this->Num_Diff(_n_xi_Min,_I_mob,_conc_Mob,delta_xi,std::bind(residual_xi_Min_tilde, _2, _3, _4, _conc_NonMin_bar, _conc_Min_bar, _vec_XiMinTilde));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde, _I_mob, _n_xi_Min_tilde, _I_NMin_bar)
//	this->Num_Diff(_n_xi_Min,_I_NMin_bar,_conc_NonMin_bar,delta_xi,std::bind(residual_xi_Min_tilde, _1, _3, _4, _conc_Mob, _conc_Min_bar, _vec_XiMinTilde));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde, _I_mob + _I_NMin_bar, _n_xi_Min_tilde, _I_min)
//	this->Num_Diff(_n_xi_Min,_n_xi_Min,_conc_Min_bar,delta_xi,std::bind(residual_xi_Min_tilde, _1, _2, _4, _conc_Mob, _conc_NonMin_bar, _vec_XiMinTilde));
//
//	// Eq. 3.60
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Kin, _I_mob)
//	this->Num_Diff(_n_xi_Kin, _I_mob, _conc_Mob, delta_xi, std::bind(residual_xi_Kin, _2, _vec_Xikin));
//
//	//Eq. 3.61
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin , 0, _n_xi_Sorp, _I_mob)
//	this->Num_Diff(_n_xi_Sorp, _I_mob, _conc_Mob, delta_xi, std::bind(residual_conc_Sorp, _2, _conc_NonMin_bar));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin , _I_mob, _n_xi_Sorp, _I_NMin_bar)
//	this->Num_Diff(_n_xi_Sorp, _I_NMin_bar, _conc_NonMin_bar, delta_xi, std::bind(residual_conc_Sorp, _1, _conc_Mob));
//
//	//Eq. 3.62
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp, 0, _n_xi_Min, _I_mob)
//	this->Num_Diff(_n_xi_Min,_I_mob, _conc_Mob, delta_xi, std::bind(residual_conc_Min, _2, _AI));
//
//	//Eq. 3.63
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min, _I_mob, _n_eta_bar, _I_NMin_bar)
//	this->Num_Diff(_n_eta_bar, _I_NMin_bar, _conc_NonMin_bar, delta_xi, std::bind(residual_Eta_bar, _2, _3, _conc_Min_bar, _vec_etabar));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min, _I_mob + _I_NMin_bar, _n_eta_bar, _I_min)
//	this->Num_Diff(_n_eta_bar, _I_min, _conc_Min_bar, delta_xi, std::bind(residual_Eta_bar, _1, _3, _conc_NonMin_bar, _vec_etabar));
//
//	//Eq. 3.64
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _I_mob , _n_xi_Kin_bar, _I_NMin_bar)
//	this->Num_Diff(_n_xi_Kin_bar, _I_NMin_bar, _conc_NonMin_bar, delta_xi, std::bind(residual_xi_KinBar_Eq, _2, _3, _conc_Min_bar, _vec_XiBarKin));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _I_mob + _I_NMin_bar, _n_xi_Kin_bar, _I_min)
//	this->Num_Diff(_n_xi_Kin_bar, _I_min, _conc_Min_bar, delta_xi, std::bind(residual_xi_KinBar_Eq, _1, _3, _conc_NonMin_bar, _vec_XiBarKin));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _I_mob + _I_NMin_bar + _I_min, _n_xi_Kin_bar, _n_xi_Kin_bar)
//	this->Num_Diff(_n_xi_Kin_bar, _n_xi_Kin_bar, _vec_XiBarKin, delta_xi, std::bind(residual_xi_KinBar_Eq, _1, _2, _conc_NonMin_bar, _conc_Min_bar));
//
//	//Eq. 3.65
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, 0, _n_xi_Kin_bar, _I_mob)
//	this->Num_Diff(_n_xi_Kin_bar, _I_mob, _conc_Mob, delta_xi, std::bind(residual_xi_KinBar_Kin, _2, _3, _4, _5, _6, _conc_NonMin_bar, _conc_Min_bar, _vec_XiBarKin, _vec_XiBarKin_old, deltaT));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _I_mob, _n_xi_Kin_bar, _I_NMin_bar)
//	this->Num_Diff(_n_xi_Kin_bar, _I_mob, _conc_NonMin_bar, delta_xi, std::bind(residual_xi_KinBar_Kin, _1, _3, _4, _5, _6, _conc_Mob, _conc_Min_bar, _vec_XiBarKin, _vec_XiBarKin_old, deltaT));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _I_mob + _I_NMin_bar, _n_xi_Kin_bar, _I_min)
//	this->Num_Diff(_n_xi_Kin_bar, _I_mob, _conc_Min_bar, delta_xi, std::bind(residual_xi_KinBar_Kin, _1, _2, _4, _5, _6, _conc_Mob, _conc_NonMin_bar, _vec_XiBarKin, _vec_XiBarKin_old, deltaT));
//
//    _mat_Jacobi.block(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _I_mob + _I_NMin_bar + _I_min, _n_xi_Kin_bar, _n_xi_Kin_bar)
//	this->Num_Diff(_n_xi_Kin_bar, _I_mob, _conc_Min_bar, delta_xi, std::bind(residual_xi_KinBar_Kin, _1, _2, _3, _5, _6, _conc_Mob, _conc_NonMin_bar, _conc_Min_bar, _vec_XiBarKin_old, deltaT));
//
//	return _mat_Jacobi;

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Jacobi Matrix: \n";
	std::cout << _mat_Jacobian << std::endl;
	// end of debugging-------------------
#endif
}


//should be in the main function to avoid repetition
//void LocalProblem::read_logK(std::vector<ogsChem::chemReactionEq*>      & list_eq_reactions)
//{
//    size_t i;
//    _vec_lnK = ogsChem::LocalVector::Zero(_J);
//    for ( i=0 ; i < list_eq_reactions.size(); i++ )
//        _vec_lnK(i) = list_eq_reactions[i]->get_ln_K();
//}


void LocalProblem::Solv_Minimization(size_t      & idx_node,
                              ogsChem::LocalMatrix & J,
                              ogsChem::LocalVector & res,
                              ogsChem::LocalVector & delta_x)
{
    size_t n_J_rows, r; 
    size_t n_R_cols;
    size_t n_V; 
    ogsChem::LocalMatrix Q, R, P; 
    ogsChem::LocalMatrix Q2, R2;
    ogsChem::LocalMatrix B; 
    ogsChem::LocalMatrix V, Vsize; 
    ogsChem::LocalVector b, b2, z; 
    ogsChem::LocalVector y1, y2, y; 
    
    b = -1.0 * res; 
    n_J_rows = J.rows(); 
    
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr_decomp; 

    // perform the QR decompostion
    qr_decomp.compute(J); 

    Q = qr_decomp.matrixQ(); 

    // TO check, is this right?
    R = qr_decomp.matrixQR(); 
    // TO check, is this right?
    P = qr_decomp.colsPermutation(); 

    // rank revealing
    r = qr_decomp.rank(); 
    
    if ( r == n_J_rows )
    {
        // using the standard direct solver
        delta_x = J.fullPivHouseholderQr().solve( b ); 
    }  // end of if
    else if ( r < n_J_rows )
    {
        // n_R_cols = size( R, 2); 
        n_R_cols = R.cols(); 
        
        // R2 = R(1:r, 1:r);
        R2 = R.topLeftCorner(r, r);
        
        // B  = R(1:r, r+1:n_R_cols);
        B = R.topRightCorner(r, n_R_cols - r);
        
        // V = R2 \ B; 
        V = R2.fullPivHouseholderQr().solve( B );
        
        // Vsize = V'*V;
        Vsize = V.transpose() * V; 
        
        // n_V = size(Vsize,1);
        n_V = Vsize.rows(); 
        
        // Q2 = Q.topLeftCorner(1:r, 1:r);
        Q2 = Q.topLeftCorner(r, r);

        // b2 = b(1:r, :);
        b2 = b.head(r); 
        
        // z = R2 \ (Q2'*b2);
        z = R2.fullPivHouseholderQr().solve( Q2.transpose() * b2 );

        // y2 = (eye(n_V)+V'*V)\(V' * z);
        ogsChem::LocalMatrix eye = ogsChem::LocalMatrix::Identity( n_V, n_V );
        y2 = ( eye + V.transpose() * V ).fullPivHouseholderQr().solve( V.transpose()*z ); 
    
        // y1 = z - V*y2; 
        y1 = z - V * y2; 

        // y = [y1;y2];
        y = ogsChem::LocalVector::Zero( y1.rows() + y2.rows() );
        y.head( y1.rows() ) = y1;
        y.tail( y2.rows() ) = y2;
        // delta_x = P* y;
        delta_x = P * y; 
    }  // end of else if

}

void LocalProblem::increment_unknown(ogsChem::LocalVector & x_old,
                                       ogsChem::LocalVector & delta_x,
                                       ogsChem::LocalVector & x_new)
{
    size_t i, _n_unknowns;
    double damp_factor, tmp_value;

    _n_unknowns = x_old.rows();
    // increment with a damping factor for everyone
    for (i=0; i<_n_unknowns; i++)
    {
    	tmp_value = -1.33*delta_x(i) / x_old(i);
        damp_factor = 1.0 / std::max(1.0, tmp_value );
        x_new(i) = x_old(i) + damp_factor * delta_x(i);
    }  // end of for

}  // end of func increment_unknown

void LocalProblem::update_minerals_conc_AI(ogsChem::LocalVector & vec_unknowns,
									      ogsChem::LocalVector & vec_AI)
{
    size_t i, idx; 
    double  cbarmin_update, phi(0.0);
    ogsChem::LocalVector conc_Mob, conc_NonMin_bar, conc_Min_bar, logConc_Mob;
    ogsChem::LocalMatrix mat_S1min_transposed;
    logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);
    mat_S1min_transposed = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(),_mat_S1min.rows());
    mat_S1min_transposed = _mat_S1min.transpose();


    // take the first section which is basis concentration
    conc_Mob    = vec_unknowns.head( _I_mob   );
    // take the non mineral or sorbed concentrations
    conc_NonMin_bar    = vec_unknowns.segment( _I_mob,_I_NMin_bar );
    //take the mineral parts
    conc_Min_bar  = vec_unknowns.tail( _I_min );


    // take the log of mobile basis concentrations
    this->cal_log_conc_vec(conc_Mob, logConc_Mob);


    // primary assessment: _AI = 1 if mineral is present and _AI = 0 if mineral is not present.
    for ( i=0; i<_I_min; i++ )
    {
    	if (conc_Min_bar(i) > 0.0)
    	{
    		vec_AI(i) = 1;
    	}
    }

    // secondary assessment: recalculate the mineral concentration and re-evaluate the _AI values.
    for ( i=0; i < _I_min; i++ )
    {
        idx = _I_mob + _I_NMin_bar + i;
        if ( vec_AI(i) == 1 )
        {
        	phi  = -_logk_min(i) + mat_S1min_transposed.row(i) * logConc_Mob;

        	if ( phi < conc_Min_bar(i))
        	{
        		// mineral is presenet
        		vec_AI(i) = 1;
        		// update mineral concentration
        		cbarmin_update = cal_cbarmin_by_constrain(i, conc_Mob, conc_NonMin_bar, conc_Min_bar);
        		vec_unknowns(idx) = cbarmin_update;
        	}  // end of if

        	if (phi >= conc_Min_bar(i))
        	{
        		vec_AI(i) = 0;
        		cbarmin_update = 0.0;
        		vec_unknowns(idx) = cbarmin_update;
        	} // end of if
        	else
            {
        			if(vec_AI(i) == 0)
        			{
        				cbarmin_update = cal_cbarmin_by_constrain(i, conc_Mob, conc_NonMin_bar, conc_Min_bar);
        				vec_unknowns(idx) = cbarmin_update;
        			} // end of if

        			vec_AI(i) = 1;

        	} // end of else

        }  // end of if AI(i)

        //vec_unknowns(idx) = cbarmin(i);
    }  // end of for loop

}  // end of function update_minerals


double LocalProblem::cal_cbarmin_by_constrain(size_t        idx_min,
                                                 ogsChem::LocalVector & conc_Mob,
                                                 ogsChem::LocalVector & conc_NonMin_bar,
                                                 ogsChem::LocalVector & conc_Min_bar)
{

    double cbarmin (0.0);
    ogsChem::LocalVector temp, conc_bar;
    ogsChem::LocalMatrix xi_sorp_bar_Ald;
    ogsChem::LocalVector xi_mobile;
    temp 				= ogsChem::LocalVector::Zero(_n_xi_Kin_bar + _n_xi_Sorp_bar + _n_xi_Min_bar);
    conc_bar 			= ogsChem::LocalVector::Zero(_n_xi_Kin_bar + _n_xi_Sorp_bar + _n_xi_Min_bar);
    xi_mobile           = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
    xi_sorp_bar_Ald     = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

    conc_bar.head(_n_xi_Kin_bar + _n_xi_Sorp_bar)  = conc_NonMin_bar;
    conc_bar.tail(_n_xi_Min_bar)   = conc_Min_bar;

    xi_mobile = _mat_c_mob_2_xi_mob * conc_Mob;

    // xi_sorp_bar
    temp = _mat_c_immob_2_xi_immob * conc_bar;
    xi_sorp_bar_Ald = _mat_Ald * temp.segment(_n_xi_Sorp_bar_li,_n_xi_Sorp_bar_ld);
    
    cbarmin = - _vec_XiMinTilde(idx_min) + xi_mobile(_n_xi_Mob + _n_xi_Sorp_bar_li + idx_min) - xi_sorp_bar_Ald(idx_min);


return cbarmin;
}  // end of function cal_cbarmin_by_total_mass

// take log of concentraiton vectors
void LocalProblem::cal_log_conc_vec(ogsChem::LocalVector & conc_Mob,
									ogsChem::LocalVector & logConc_Mob)
{
	double tmp_x;
	std::size_t i;


	for (i = 0; i < _I_mob; i++)
	{
		tmp_x    = conc_Mob(i);
		logConc_Mob(i)  = std::log(tmp_x);
	}
#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "conc_Mob Vector: \n";
	std::cout << conc_Mob << std::endl;

	std::cout << "logConc_Mob Vector: \n";
	std::cout << logConc_Mob << std::endl;

	// end of debugging-------------------
#endif
}

// Eq. 3.55
void LocalProblem::residual_conc_Mob(ogsChem::LocalVector & conc_Mob,
									 ogsChem::LocalVector & vec_residual)
{
	 ogsChem::LocalVector logConc_Mob;
	 logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);


	 this->cal_log_conc_vec(conc_Mob, logConc_Mob);
	vec_residual.head(_n_xi_Mob) 	= - _logk_mob + _mat_S1mob.transpose() * logConc_Mob;

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "_logk_mob Vector: \n";
	std::cout << _logk_mob << std::endl;

	std::cout << "_mat_S1mob: \n";
	std::cout << _mat_S1mob << std::endl;

	std::cout << "logConc_Mob vector: \n";
	std::cout << logConc_Mob << std::endl;

	// end of debugging-------------------
#endif

}

// Eq. 3.56
void LocalProblem::residual_Eta(ogsChem::LocalVector & conc_Mob,
		                        ogsChem::LocalVector & vec_residual)
{
	// eta acts as a constrain
	vec_residual.segment(_n_xi_Mob, _n_eta)  = - _vec_eta + _mat_c_mob_2_eta_mob * conc_Mob;
}

// Eq. 3.57 - 58
void LocalProblem::residual_xi_Sorp_tilde(ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_NonMin_bar,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector                 vec_XiSorp;
	ogsChem::LocalVector                 vec_XiSorpBarLI;
	ogsChem::LocalVector                 conc_bar;
	conc_bar = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);

	conc_bar.head(_n_xi_Sorp_bar) = conc_NonMin_bar;
	conc_bar.tail(_n_xi_Min) 	  = conc_Min_bar;

	vec_XiSorp        = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_XiSorpBarLI   = _mat_c_immob_2_xi_immob * conc_bar;
	vec_residual.segment(_n_xi_Mob + _n_eta, _n_xi_Sorp_tilde) = - _vec_XiSorpTilde + vec_XiSorp.segment(_n_xi_Mob,_n_xi_Sorp_bar_li ) - vec_XiSorpBarLI.segment(0,_n_xi_Sorp_bar_li);
}

// Eq. 3.59
void LocalProblem::residual_xi_Min_tilde(ogsChem::LocalVector & conc_Mob,
										 ogsChem::LocalVector & conc_NonMin_bar,
										 ogsChem::LocalVector & conc_Min_bar,
										 ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector                 vec_XiMin;
	ogsChem::LocalVector                 vec_XiMinBar;
	ogsChem::LocalVector                 conc_bar;
	conc_bar 		= ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);
	vec_XiMin 		= ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	vec_XiMinBar 	= ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);

	conc_bar.head(_n_xi_Sorp_bar)    = conc_NonMin_bar;
	conc_bar.tail(_n_xi_Min)         = conc_Min_bar;

	vec_XiMin        = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_XiMinBar     = _mat_c_immob_2_xi_immob * conc_bar;
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde, _n_xi_Min_tilde) = - _vec_XiMinTilde + vec_XiMin.segment(_n_xi_Mob + _n_xi_Sorp_bar_li, _n_xi_Min) - conc_Min_bar - (_mat_Ald * vec_XiMinBar.tail(_n_xi_Sorp_bar_ld));
}

// Eq. 3.60
void LocalProblem::residual_xi_Kin(ogsChem::LocalVector & conc_Mob,
								   ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector   conc_tmp;
	conc_tmp 		= ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);

	conc_tmp     = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde , _n_xi_Kin)   = - _vec_Xikin + conc_tmp.segment(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin);
}

// Eq. 3.61
void LocalProblem::residual_conc_Sorp(ogsChem::LocalVector & conc_Mob,
									  ogsChem::LocalVector & conc_NonMin_bar,
									  ogsChem::LocalVector & vec_residual)
{

	ogsChem::LocalVector   conc_tmp;
	ogsChem::LocalVector logConc_Mob;
	logConc_Mob  = ogsChem::LocalVector::Zero(_I_mob + _I_NMin_bar);
	conc_tmp	 = ogsChem::LocalVector::Zero(_I_mob + _I_NMin_bar);


	//TODO initialize log k values
	conc_tmp.head(_I_mob)  = conc_Mob ;
	conc_tmp.tail(_I_NMin_bar) = conc_NonMin_bar;

	this->cal_log_conc_vec(conc_tmp, logConc_Mob);
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +_n_xi_Kin, _n_xi_Sorp)  = - _logk_sorp + _mat_Ssorp.transpose() * logConc_Mob;
}

// Eq. 3.62
void LocalProblem::residual_conc_Min(ogsChem::LocalVector & conc_Mob,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & vec_AI,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & vec_residual)
{
	size_t        i, idx;
	ogsChem::LocalMatrix   mat_S1minT;
	ogsChem::LocalVector logConc_Mob;
	logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);
	mat_S1minT 	= ogsChem::LocalMatrix::Zero(_mat_S1min.cols(), _mat_S1min.rows());

	idx  = _n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp;

	mat_S1minT = _mat_S1min.transpose();

	for (i=0; i < _n_xi_Min; i++)
	{
		if (vec_AI(i) == 1)
		{
		this->cal_log_conc_vec(conc_Mob, logConc_Mob);
		vec_residual(idx + i) 	= - _logk_min(i) + mat_S1minT.row(i) * logConc_Mob;
		}

	}
}

// Eq. 3.63
void LocalProblem::residual_Eta_bar(ogsChem::LocalVector & conc_NonMin_bar,
									ogsChem::LocalVector & conc_Min_bar,
									ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector      conc_bar;
	conc_bar  = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);

	conc_bar.head(_I_NMin_bar)    = conc_NonMin_bar;
	conc_bar.tail(_I_min)         = conc_Min_bar;

	// etabar acts as a constrain
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min, _n_eta_bar) 	  = - _vec_etabar + (_mat_c_immob_2_eta_immob * conc_bar);
}

// Eq. 3.64
void LocalProblem::residual_xi_KinBar_Eq(ogsChem::LocalVector & conc_NonMin_bar,
										 ogsChem::LocalVector & conc_Min_bar,
										 ogsChem::LocalVector & Xi_Kin_bar,
										 ogsChem::LocalVector & vec_residual)
{

	ogsChem::LocalVector   conc_tmp;
	ogsChem::LocalVector   conc_bar;
	conc_bar  = ogsChem::LocalVector::Zero(_I_NMin_bar + _n_xi_Min);
	conc_tmp  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);

	conc_bar.head(_I_NMin_bar)    = conc_NonMin_bar;
	conc_bar.tail(_n_xi_Min)         = conc_Min_bar;

	conc_tmp     = _mat_c_immob_2_xi_immob * conc_bar;
	//XiBarKin is the total mass constrain
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _n_xi_Kin_bar)  = - Xi_Kin_bar + conc_tmp.segment(_n_xi_Sorp_bar + _n_xi_Min_bar, _n_xi_Kin_bar);
}

// Eq. 3.65
void LocalProblem::residual_xi_KinBar_Kin(ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_NonMin_bar,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & Xi_Kin_bar,
										  ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector   conc, vec_rateKin;
	conc 		 = ogsChem::LocalVector::Zero(_n_Comp);
	vec_rateKin  = ogsChem::LocalVector::Zero(_n_xi_Kin);

	//TODO get theta_waterContent
	double theta_waterContent (0.3);

	conc.head	(_I_mob)  					 = conc_Mob;
	conc.segment(_I_mob, _I_NMin_bar) 		 = conc_NonMin_bar;
	conc.tail	(_n_xi_Min_bar) 			 = conc_Min_bar;

	this->reaction_rates(conc, vec_rateKin);

	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _n_xi_Kin_bar)
			      = ((theta_waterContent * Xi_Kin_bar - (theta_waterContent * _vec_XiBarKin_old)) / deltaT) - (theta_waterContent * _mat_A2kin * vec_rateKin);
}

//problem specific reaction rates
void LocalProblem::reaction_rates(ogsChem::LocalVector & conc,
								  ogsChem::LocalVector & vec_rateKin)
{
//sub_A  = conc.segment(0,0);
//sub_B  = conc.segment(1,0);
//sub_C  = conc.segment(2,0);
//biomass      = conc.segment(3,0);
//
//// max reaction rate
//umax = umax;
//
//// monod coefficients
//k_subA = k_subA;
//k_subB = k_subB;
//
//// decay rates
//kdec = kdec;
//
//R.head(0,0) = umax * sub_A / (k_subA + sub_A) * sub_B /(k_subB + sub_B) * biomass;
//R.tail (1,0) =  R(0,1) - kdec * biomass;
}

//template <typename T>
//ogsChem::LocalMatrix LocalProblem::Num_Diff(size_t & Row,
//							std::size_t & Col,
//							ogsChem::LocalVector & _vec_unknown,
//							double & delta_xi,
//							std::function<ogsChem::LocalVector(T))> f)
//{
//    size_t i;
//    res = ogsChem::LocalMatrix::Zero(Row, Col);
//	for (i = 0; i< Col; i++)
//	{
//		x  = _vec_unknown;
//		x(i) = x(i) +  delta_xi;
//		res.Col(i) = ( f(x) - f( _vec_unknown)) / delta_xi;
//	}
//
//	return res;
////	return f(_vec_unknown);
//
//}
