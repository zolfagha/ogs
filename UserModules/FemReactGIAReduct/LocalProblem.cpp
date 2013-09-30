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
//#define _DEBUG

LocalProblem::LocalProblem(ogsChem::chemReductionGIA* ReductionGIA, MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>* sbs, Local_ODE_Xi_immob_GIA* local_ode_xi_immob_GIA)
	: _n_Comp(ReductionGIA->get_n_Comp()), _ReductionGIA(ReductionGIA), _sbs(sbs), _local_ode_xi_immob_GIA(local_ode_xi_immob_GIA), _I_mob(_ReductionGIA->get_n_Comp_mob()), _I_sorp(_ReductionGIA->get_n_Comp_sorb()), _I_min(_ReductionGIA->get_n_Comp_min()), _n_xi_Kin_bar(_ReductionGIA->get_n_xi_Kin_bar())
	, _n_xi_Mob(_ReductionGIA->get_n_xi_Mob()), _n_eta(_ReductionGIA->get_n_eta()), _n_eta_bar(_ReductionGIA->get_n_eta_bar()), _n_xi_Sorp_tilde(_ReductionGIA->get_n_xi_Sorp_tilde()), _n_xi_Min_tilde(_ReductionGIA->get_n_xi_Min_tilde())
    , _n_xi_Sorp(_ReductionGIA->get_n_xi_Sorp()), _n_xi_Min(_ReductionGIA->get_n_xi_Min()), _n_xi_Sorp_bar_li(_ReductionGIA->get_n_xi_Sorp_bar_li()), _n_xi_Sorp_bar_ld(_ReductionGIA->get_n_xi_Sorp_bar_ld()), _n_xi_Kin(_ReductionGIA->get_n_xi_Kin())
	, _mat_c_mob_2_xi_mob(_ReductionGIA->get_matrix_C2Xi()), _mat_c_immob_2_xi_immob(_ReductionGIA->get_matrix_Cbar2XiBar()), _n_xi_Sorp_bar(_ReductionGIA->get_n_xi_Sorp_bar()), _mat_S1min (_ReductionGIA->get_matrix_S1min())
    , _mat_Ald(_ReductionGIA->get_matrix_Ald()),_mat_Ssorp(_ReductionGIA->get_matrix_Ssorp()), _mat_A2kin(_ReductionGIA->get_matrix_A2kin()), _mat_S1mob(_ReductionGIA->get_matrix_S1mob()),_mat_c_mob_2_eta_mob(_ReductionGIA->get_matrix_C2Eta())
    , _logk_mob(_ReductionGIA->get_logk_mob()), _logk_sorp(_ReductionGIA->get_logk_sorp()), _logk_min(_ReductionGIA->get_logk_min()), _n_xi_Min_bar(_ReductionGIA->get_n_xi_Min_bar()), _I_NMin_bar(_ReductionGIA->get_n_Comp_NMin_bar())
    , _mat_c_immob_2_eta_immob(_ReductionGIA->get_matrix_C2EtaBar()), _list_kin_reactions(_ReductionGIA->get_list_kin_reactions()), _J_tot_kin(_ReductionGIA->get_n_xi_Kin_total()), _n_xi_global(_ReductionGIA->get_n_xi_global())
	, _n_xi_local(_ReductionGIA->get_n_xi_local())

{

}

//LocalProblem::~LocalProblem()
//{}

void LocalProblem::solve_LocalProblem_Newton_LineSearch(std::size_t & node_idx,
														double dt,
														const double iter_tol,
														const double rel_tol,
														const double max_iter, 
														ogsChem::LocalVector & x,
														ogsChem::LocalVector & vec_eta, 
														ogsChem::LocalVector & vec_etabar,
														ogsChem::LocalVector & vec_XiSorpTilde, 
														ogsChem::LocalVector & vec_XiMinTilde, 
														ogsChem::LocalVector & vec_Xikin, 
														ogsChem::LocalVector & vec_XiBarKin, 
														ogsChem::LocalVector & vec_XiBarKin_old)
{
    ogsChem::LocalVector x_new, vec_residual, vec_AI;
    ogsChem::LocalVector dx;

    // TODO initialize water content and deltaT

    // number of iterations
    size_t j, iter;
    const double alpha (0.5);
    double d_norm(0.0), d1_norm(0.0);

	// HS: notice that the number of unknowns equals to _n_Comps
    x_new          = ogsChem::LocalVector::Zero( _n_Comp );
	dx             = ogsChem::LocalVector::Ones( _n_Comp );
	_mat_Jacobian  = ogsChem::LocalMatrix::Zero( _n_Comp, _n_Comp );
	vec_residual   = ogsChem::LocalVector::Zero( _n_Comp );
    // initialize the _AI vector
    vec_AI		   = ogsChem::LocalVector::Zero( _I_min );

	// vec_tot_mass_constrain contains xi global and eta mobile which acts as a total mass constrain for the local problem.
	_vec_eta            = vec_eta;
	_vec_etabar         = vec_etabar;
	_vec_XiSorpTilde    = vec_XiSorpTilde;
	_vec_XiMinTilde     = vec_XiMinTilde;
	_vec_Xikin          = vec_Xikin;
	_vec_XiBarKin       = vec_XiBarKin;
	_vec_XiBarKin_old   = vec_XiBarKin_old;

#ifdef _DEBUG
    // // debugging--------------------------
//     std::cout << "======================================== \n";
//     std::cout << "vec_tot_mass_constrain: \n";
//     std::cout << vec_tot_mass_constrain << std::endl;
//     std::cout << "x: \n";
//     std::cout << x << std::endl;
//     std::cout << "======================================== \n";
    // // end of debugging-------------------
#endif

    // start solving the system
    iter = 0; 

    // now updating the saturation index and minerals
    if(_n_xi_Min != 0)
    this->update_minerals_conc_AI( x, vec_AI );

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "x Vector: \n";
//	std::cout << x << std::endl;
//
//	std::cout << "vec_AI Vector: \n";
//	std::cout << vec_AI << std::endl;
	// end of debugging-------------------
#endif


	// update the value of xikinbar in the vector of unknowns(x_new)
    //this->ODE_solver(dt, x);
    // evaluate the residual
	this->calc_residual(dt, x, _vec_XiBarKin_old, vec_residual, _vec_XiBarKin);

    //save the update value of xi kin bar
    //vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_eta_bar, _n_xi_Kin_bar) = _vec_XiBarKin;
    // evaluate norm of residual vector
    d_norm = vec_residual.norm();

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Residual Vector: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

	// save the previous values
	x_new  =  x;

    while (true)
    {
        #ifdef _DEBUG
            // display the residual
            // std::cout << "Iteration #" << iter << "||res|| = " << d_norm << "||delta_x|| = " << dx.norm() << std::endl;
        #endif
        // convergence criteria
        if ( d_norm < iter_tol )
        {
            #ifdef _DEBUG
 //                std::cout << "Newton iteration successfully converged!\n";
            #endif

            break;  // break the loop
        }
        else if ( dx.norm() < rel_tol )
        {
            #ifdef _DEBUG
 //           std::cout << "Warning, Newton iteration stagnent on Node #" << node_idx << "! Exit the iteration!\n" ;
            #endif

            break;  // break the loop
        }
        else if ( iter > max_iter )
        {
            #ifdef _DEBUG
//            std::cout << "ERROR! Node #" << node_idx  << "Newton iterationan does not converge! Simulation stops!\n";
            #endif

            return; // stop the program
        }
        // form Jacobian matrix
		this->calc_Jacobian(dt, x, vec_residual, _vec_XiBarKin);

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "_mat_Jacobian: \n";
	std::cout << _mat_Jacobian << std::endl;
	// end of debugging-------------------
#endif

        // solving for increment
        this->Solv_Minimization( node_idx, _mat_Jacobian, vec_residual, dx );

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "dx Vector: \n";
//	std::cout << dx << std::endl;
	// end of debugging-------------------
#endif


        // increment of unknowns
        this->increment_unknown( x, dx, x_new ); 

#ifdef _DEBUG
	// debugging--------------------------
//	 std::cout << "x_new Vector: \n";
//	 std::cout << x_new << std::endl;
	// end of debugging-------------------
#endif

    	// updating the saturation index and minerals
        if(_n_xi_Min != 0)
		this->update_minerals_conc_AI( x_new, vec_AI );

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "x_new: \n";
//	std::cout << x_new << std::endl;
	// end of debugging-------------------
#endif
		// update the value of xikinbar in the vector of unknowns(x_new)
		this->ODE_solver(dt, x_new, _vec_XiBarKin);
        // evaluate residual with x_new
		this->calc_residual(dt, x_new, _vec_XiBarKin_old, vec_residual, _vec_XiBarKin);


#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "x_new Vector: \n";
	std::cout << x_new << std::endl;

	std::cout << "vec_residual: \n";
	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

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
            if(_n_xi_Min != 0)
            this->update_minerals_conc_AI( x_new, vec_AI );
            // update the value of xikinbar in the vector of unknowns(x_new)
			this->ODE_solver(dt, x_new, _vec_XiBarKin);
			// evaluate residual with x_new
			this->calc_residual(dt, x_new, _vec_XiBarKin_old, vec_residual, _vec_XiBarKin);

		#ifdef _DEBUG
        	std::cout << "vec_residual: \n";
        	std::cout << vec_residual << std::endl;
            // display the residual
            // std::cout << "Line Search Iteration #" << iter << "||res|| = " << d1_norm << "||delta_x|| = " << dx.norm() << std::endl;
		#endif

            j++;
        }  // end of while
        d_norm = d1_norm; 
        x = x_new; 

        // increase the iteration count
        iter++; 
    }
}

void LocalProblem::ODE_solver(double dt,
							  ogsChem::LocalVector & vec_unknowns,
							  ogsChem::LocalVector & vec_Xi_Kin_bar)
{
	ogsChem::LocalVector conc_Mob 		    = ogsChem::LocalVector::Zero(_I_mob);
	ogsChem::LocalVector ln_conc_Mob	    = ogsChem::LocalVector::Zero(_I_mob);
	ogsChem::LocalVector conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
	ogsChem::LocalVector ln_conc_NonMin_bar = ogsChem::LocalVector::Zero(_I_NMin_bar);
	ogsChem::LocalVector vec_xi_kin_rate    = ogsChem::LocalVector::Zero(_J_tot_kin);
	ogsChem::LocalVector conc_Min_bar	    = ogsChem::LocalVector::Zero(_I_min);

    ln_conc_Mob         = vec_unknowns.segment(0, _I_mob );
    //NonMin refers to sorbed and kinetic concentrations.
    ln_conc_NonMin_bar  = vec_unknowns.segment(_I_mob, _I_NMin_bar );
    conc_Min_bar        = vec_unknowns.segment(_I_mob + _I_NMin_bar, _I_min ); // nonmineral conc = xi kin bar only in this case

    MathLib::LocalVector  loc_eta        	= MathLib::LocalVector::Zero( _n_eta );
    MathLib::LocalVector  loc_eta_bar     	= MathLib::LocalVector::Zero( _n_eta_bar );
    MathLib::LocalVector  loc_xi_global     = MathLib::LocalVector::Zero( _n_xi_global );
    MathLib::LocalVector  loc_xi_local      = MathLib::LocalVector::Zero( _n_xi_local );
    MathLib::LocalVector  local_conc        = MathLib::LocalVector::Zero( _n_Comp );
    MathLib::LocalVector  vec_conc          = MathLib::LocalVector::Zero( _n_Comp );



	vec_conc   	      =  vec_unknowns.head (_n_Comp);

    // convert the ln mobile conc to mobile conc
    ln_conc_Mob 	  =  vec_conc.head(_I_mob);
    this->cal_exp_conc_vec(_I_mob,ln_conc_Mob, conc_Mob);

    conc_Min_bar	 =  vec_conc.tail(_I_min);
     // convert the log nonmineral conc to nonmineral conc
    ln_conc_NonMin_bar  = vec_conc.segment(_I_mob, _I_NMin_bar);
    this->cal_exp_conc_vec(_I_NMin_bar, ln_conc_NonMin_bar, conc_NonMin_bar);

    local_conc.head(_I_mob)				 	  =  conc_Mob;
    local_conc.segment(_I_mob, _I_NMin_bar)   =  conc_NonMin_bar;
    local_conc.tail(_I_min) 			      =  conc_Min_bar;

    _ReductionGIA->Conc2EtaXi(local_conc, loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local);
    //using ode solver for vec_XiBarKin
	// get the right reference values to ODE RHS function
//	this->_local_ode_xi_immob_GIA->update_eta_xi(loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local, Xi_Kin_bar);
//	vec_xi_kin_rate = (*_local_ode_xi_immob_GIA)(dt, Xi_Kin_bar);

	this->_local_ode_xi_immob_GIA->update_eta_xi(loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local);
	vec_xi_kin_rate = (*_local_ode_xi_immob_GIA)(dt, vec_Xi_Kin_bar);

    //_sbs->set_y(Xi_Kin_bar);
	_sbs->set_y(conc_NonMin_bar);
    _sbs->set_dydx(vec_xi_kin_rate);

	// solve the local ODE problem for xi kin bar
    _sbs->step(dt, _local_ode_xi_immob_GIA);
    //vec_Xi_Kin_bar_new = _sbs->get_y();
    vec_Xi_Kin_bar = _sbs->get_y();


}  // end of function ode solver


void LocalProblem::calc_residual(double dt,
								 ogsChem::LocalVector & vec_unknowns,
								 ogsChem::LocalVector & vec_xi_Kin_bar_old, 
								 ogsChem::LocalVector & vec_residual,
								 ogsChem::LocalVector & vec_Xi_Kin_bar)
{
	ogsChem::LocalVector ln_conc_Mob, ln_conc_NonMin_bar, conc_Mob, conc_NonMin_bar, conc_Min_bar;
	conc_Mob 		    = ogsChem::LocalVector::Zero(_I_mob);
	conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);

	ln_conc_Mob = vec_unknowns.segment(0, _I_mob);
    //NonMin refers to sorbed and kinetic concentrations.
	ln_conc_NonMin_bar = vec_unknowns.segment(_I_mob, _I_NMin_bar);
	conc_Min_bar = vec_unknowns.segment(_I_mob + _I_NMin_bar, _I_min);

    // convert the ln mobile conc to mobile conc
     this->cal_exp_conc_vec(_I_mob, ln_conc_Mob, conc_Mob);
     // convert the ln immobile conc to immobile conc
     this->cal_exp_conc_vec(_I_NMin_bar, ln_conc_NonMin_bar, conc_NonMin_bar);

    // Eq. 3.55
    this->residual_conc_Mob			(ln_conc_Mob, vec_residual);
    // Eq. 3.56
    this->residual_Eta				(conc_Mob, vec_residual);

    if(_n_xi_Min != 0)
    {
    // Eq. 3.59
    this->residual_xi_Min_tilde		(conc_Mob, conc_NonMin_bar, conc_Min_bar, vec_residual);
    // Eq. 3.62
    this->residual_conc_Min			(ln_conc_Mob, conc_Min_bar, vec_residual);
    }

    if(_n_xi_Sorp != 0)
    {
    // Eq. 3.57 - 58
    this->residual_xi_Sorp_tilde	(conc_Mob, conc_NonMin_bar, conc_Min_bar, vec_residual);
    // Eq. 3.61
    this->residual_conc_Sorp		(ln_conc_Mob, ln_conc_NonMin_bar, vec_residual);
    }
    // Eq. 3.63
    this->residual_Eta_bar			(conc_NonMin_bar, conc_Min_bar, vec_residual);

    if(_n_xi_Kin != 0)
    {
    	// Eq. 3.60
    	this->residual_xi_Kin			(conc_Mob, vec_residual);


    	// Eq. 3.64
    	this->residual_xi_KinBar_Eq		(conc_NonMin_bar, conc_Min_bar, vec_Xi_Kin_bar, vec_residual); //RZ: solve xi kin bar using ode solver after solving local problem.
    	// Eq. 3.65
    	//this->residual_xi_KinBar_Kin	(dt,conc_Mob, conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar, vec_residual);  //RZ: solve xi kin bar using ode solver after solving local problem.
    }

}  // end of function calc_residual

void LocalProblem::calc_Jacobian(double dt,
								 ogsChem::LocalVector & vec_x,
//							     ogsChem::LocalVector & vec_AI,
							     ogsChem::LocalVector & vec_residual,
							     ogsChem::LocalVector & vec_Xi_Kin_bar)
{
	//const double delta_xi = 1.0e-8;  //calcite example
	const double delta_xi = 1.0e-10;    //monod2d
    int i;
    ogsChem::LocalVector vec_x_incremented, vec_residual_incremented;
    vec_residual_incremented = vec_residual;

	for (i = 0; i < _mat_Jacobian.cols(); i++)
	{
		vec_x_incremented  = vec_x;

//		// numerical protection
//		if( vec_x_incremented.norm() < 1.0e-16)
		{
			vec_x_incremented(i) += delta_xi * vec_x_incremented.norm();
			this->calc_residual(dt, vec_x_incremented, _vec_XiBarKin_old, vec_residual_incremented, vec_Xi_Kin_bar);
			_mat_Jacobian.col(i) = (vec_residual_incremented - vec_residual ) / (delta_xi * vec_x_incremented.norm());

		}
//		else
//		{
//		vec_x_incremented(i) += delta_xi;
//
//		this->calc_residual(vec_x_incremented,vec_AI,vec_residual_incremented);
//		_mat_Jacobian.col(i) = (vec_residual_incremented - vec_residual ) / delta_xi;
//		}
	}



#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "Jacobi Matrix: \n";
	// std::cout << _mat_Jacobian << std::endl;
	// end of debugging-------------------
#endif
}


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
//    size_t i, _n_unknowns;
//    double damp_factor, tmp_value;
//    ogsChem::LocalVector ln_conc_Mob_old, ln_conc_NonMin_bar_old, conc_Min_bar_old, Xi_Kin_bar_old, delta_x_mob, delta_x_nonmin, delta_x_min, delta_x_kin
//                         , ln_conc_Mob_new, ln_conc_NonMin_bar_new, Xi_Kin_bar_new, conc_Min_bar_new;
//
//    ln_conc_Mob_old 	    = ogsChem::LocalVector::Zero(_I_mob);
//    ln_conc_NonMin_bar_old 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
//    conc_Min_bar_old 	    = ogsChem::LocalVector::Zero(_I_min);
//    Xi_Kin_bar_old 			= ogsChem::LocalVector::Zero(_n_xi_Kin_bar);
//    delta_x_mob 			= ogsChem::LocalVector::Zero(_I_mob);
//    delta_x_nonmin 			= ogsChem::LocalVector::Zero(_I_NMin_bar);
//    delta_x_min 			= ogsChem::LocalVector::Zero(_I_min);
//    delta_x_kin 			= ogsChem::LocalVector::Zero(_n_xi_Kin_bar);
//    ln_conc_Mob_new 		= ogsChem::LocalVector::Zero(_I_mob);
//    ln_conc_NonMin_bar_new 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
//    conc_Min_bar_new 		= ogsChem::LocalVector::Zero(_I_min);
//    Xi_Kin_bar_new 			= ogsChem::LocalVector::Zero(_n_xi_Kin_bar);

      x_new = x_old + delta_x;



//    ln_conc_Mob_old     	= x_old.head(_I_mob );
//    ln_conc_NonMin_bar_old  = x_old.segment(_I_mob, _I_NMin_bar );
//    conc_Min_bar_old     	= x_old.segment(_I_mob + _I_NMin_bar, _I_min );
//    Xi_Kin_bar_old       	= x_old.tail(_n_xi_Kin_bar);
//
//    delta_x_mob     	= delta_x.head(_I_mob );
//    delta_x_nonmin 	    = delta_x.segment(_I_mob, _I_NMin_bar );
//    delta_x_min     	= delta_x.segment(_I_mob + _I_NMin_bar, _I_min );
//    delta_x_kin       	= delta_x.tail(_n_xi_Kin_bar);
//
//    ln_conc_Mob_new         =  ln_conc_Mob_old        + delta_x_mob;
//    ln_conc_NonMin_bar_new  =  ln_conc_NonMin_bar_old + delta_x_nonmin;
//    Xi_Kin_bar_new          =  Xi_Kin_bar_old         + Xi_Kin_bar_old;
//    // increment with a damping factor for minerals
//    for (i=0; i<_I_min; i++)
//    {
//    	tmp_value = -1.33*delta_x_min(i) / conc_Min_bar_old(i);
//        damp_factor = 1.0 / std::max(1.0, tmp_value );
//        conc_Min_bar_new(i) = conc_Min_bar_old(i) + damp_factor * delta_x_min(i);
//    }  // end of for
//
//    x_new.head(_I_mob ) 						 = ln_conc_Mob_new;
//    x_new.segment(_I_mob, _I_NMin_bar ) 		 = ln_conc_NonMin_bar_new;
//    x_new.segment(_I_mob + _I_NMin_bar, _I_min ) = conc_Min_bar_new;
//    x_new.tail(_n_xi_Kin_bar)  				     = Xi_Kin_bar_new;

}  // end of func increment_unknown


void LocalProblem::update_minerals_conc_AI(ogsChem::LocalVector & vec_unknowns,
									       ogsChem::LocalVector & vec_AI)
{
    size_t i, idx; 
    double  phi(0.0);
    ogsChem::LocalVector ln_conc_Mob = ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalVector ln_conc_NonMin_bar = ogsChem::LocalVector::Zero(_I_NMin_bar);
    ogsChem::LocalVector conc_Min_bar = ogsChem::LocalVector::Zero(_I_min);
    //ogsChem::LocalMatrix mat_S1min_transposed;
    //logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalMatrix  mat_S1min_transposed = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(),_mat_S1min.rows());
    mat_S1min_transposed = _mat_S1min.transpose();


    // take the first section which is basis concentration
    ln_conc_Mob    = vec_unknowns.head( _I_mob   );
    // take the non mineral concentrations
    ln_conc_NonMin_bar    = vec_unknowns.segment( _I_mob,_I_NMin_bar );
    //take the mineral parts
    conc_Min_bar  = vec_unknowns.tail( _I_min );


    // _AI = 1 if mineral is present and _AI = 0 if mineral is not present.
    for ( i=0; i < _I_min; i++ )
    {
    	idx = _I_mob + _I_NMin_bar + i;

    	phi  = -_logk_min(i) + mat_S1min_transposed.row(i) * ln_conc_Mob;


    	// if mineral concentration  >= phi ; mineral is present; saturated case; precipitate the mineral.
    	if (conc_Min_bar(i) >= phi)
    	{
    		vec_AI(i) = 1;
    	}// end of if
    	// if mineral concentration < phi : mineral is NOT present; under saturated case; dissolve the mineral
    	else
    	{
    		vec_AI(i) = 0;
    		conc_Min_bar(i) = 0.0;
    	}// end of else


    	// if mineral is present, calculate mineral concentration
    	if	(vec_AI(i) == 1)
    	{
    		// update mineral concentration
    		conc_Min_bar(i) = cal_cbarmin_by_constrain(i, ln_conc_Mob, ln_conc_NonMin_bar, conc_Min_bar);
    		vec_unknowns(idx) = conc_Min_bar(i);
    	}  // end of if

    	if (phi > conc_Min_bar(i))
    	{
    		vec_AI(i) = 0;
    		conc_Min_bar(i) = 0.0;
    		vec_unknowns(idx) = conc_Min_bar(i);
    	} // end of if
    	else
    	{
    		if(vec_AI(i) == 0)
    		{
    			conc_Min_bar(i) = cal_cbarmin_by_constrain(i, ln_conc_Mob, ln_conc_NonMin_bar, conc_Min_bar);
    			vec_unknowns(idx) = conc_Min_bar(i);
    		} // end of if

    		vec_AI(i) = 1;

    	} // end of else

    }  // end of for loop

}  // end of function update_minerals


double LocalProblem::cal_cbarmin_by_constrain(size_t        idx_min,
                                                 ogsChem::LocalVector & ln_conc_Mob,
                                                 ogsChem::LocalVector & ln_conc_NonMin_bar,
                                                 ogsChem::LocalVector & conc_Min_bar)
{

    double cbarmin (0.0);
    ogsChem::LocalVector temp 				= ogsChem::LocalVector::Zero(_n_xi_Kin_bar + _n_xi_Sorp_bar + _n_xi_Min_bar);
    ogsChem::LocalVector conc_bar 			= ogsChem::LocalVector::Zero(_n_xi_Kin_bar + _n_xi_Sorp_bar + _n_xi_Min_bar);
    ogsChem::LocalVector Conc_Mob 			= ogsChem::LocalVector::Zero(_I_mob);
    ogsChem::LocalVector conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
    ogsChem::LocalVector xi_mobile           = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
    ogsChem::LocalMatrix xi_sorp_bar_Ald     = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

    // convert the ln mobile conc to mobile conc
    this->cal_exp_conc_vec(_I_mob, ln_conc_Mob, Conc_Mob);
    xi_mobile = _mat_c_mob_2_xi_mob * Conc_Mob;

    // convert the ln nonmineral immobile conc to nonmineral immobile conc
    this->cal_exp_conc_vec(_I_NMin_bar, ln_conc_NonMin_bar, conc_NonMin_bar);

    conc_bar.head(_I_NMin_bar)   = conc_NonMin_bar;
    conc_bar.tail(_n_xi_Min_bar) = conc_Min_bar;

    // xi_sorp_bar
    temp 			= _mat_c_immob_2_xi_immob * conc_bar;
    xi_sorp_bar_Ald = _mat_Ald * temp.segment(_n_xi_Sorp_bar_li,_n_xi_Sorp_bar_ld);
    
    cbarmin 		= - _vec_XiMinTilde(idx_min) + xi_mobile(_n_xi_Mob + _n_xi_Sorp_bar_li + idx_min) - xi_sorp_bar_Ald(idx_min);


return cbarmin;
}  // end of function cal_cbarmin_by_total_mass


// take log of concentration vectors
void LocalProblem::cal_ln_conc_vec(size_t                 idx_size,
								   ogsChem::LocalVector & conc_Mob,
								   ogsChem::LocalVector & ln_conc_Mob)
{
	double tmp_x;
	std::size_t i;


	for (i = 0; i < idx_size; i++)
	{
		tmp_x    = conc_Mob(i);
		ln_conc_Mob(i)  = std::log(tmp_x);
	}

}

// take exponential of natural log concentration vectors
void LocalProblem::cal_exp_conc_vec(size_t    		      idx_size,
									ogsChem::LocalVector & ln_conc,
									ogsChem::LocalVector & Conc)
{
	double tmp_x;
	std::size_t i;


	for (i = 0; i < idx_size; i++)
	{
		tmp_x    	 = ln_conc(i);
		Conc(i)  = std::exp(tmp_x);
	}

}

// Eq. 3.55
void LocalProblem::residual_conc_Mob(ogsChem::LocalVector & ln_conc_Mob,
									 ogsChem::LocalVector & vec_residual)
{

	vec_residual.head(_n_xi_Mob) 	= - _logk_mob + _mat_S1mob.transpose() * ln_conc_Mob;

}

// Eq. 3.56
void LocalProblem::residual_Eta(ogsChem::LocalVector & conc_Mob,
		                        ogsChem::LocalVector & vec_residual)
{
	// eta acts as a constrain
	vec_residual.segment(_n_xi_Mob, _n_eta)  = -1.0 * _vec_eta + _mat_c_mob_2_eta_mob * conc_Mob;

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "_vec_eta: \n";
	//std::cout << _vec_eta << std::endl;

	//std::cout << "conc_Mob: \n";
	//std::cout << conc_Mob << std::endl;

	//std::cout << "_mat_c_mob_2_eta_mob: \n";
	//std::cout << _mat_c_mob_2_eta_mob << std::endl;

	//std::cout << "vec_residual: \n";
	//std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

}

// Eq. 3.57 - 58
void LocalProblem::residual_xi_Sorp_tilde(ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_NonMin_bar,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector vec_XiSorp      = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	ogsChem::LocalVector vec_XiSorpBarLI = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector conc_bar        = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);;

	conc_bar.head(_I_NMin_bar) 	  = conc_NonMin_bar;
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
	ogsChem::LocalVector  vec_XiMin     = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	ogsChem::LocalVector  vec_XiMinBar  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector conc_bar 		= ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);
	ogsChem::LocalVector A 				= ogsChem::LocalVector::Zero(_n_xi_Min);
	ogsChem::LocalVector B 				= ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

	conc_bar.head(_I_NMin_bar)    = conc_NonMin_bar;
	conc_bar.tail(_n_xi_Min)         = conc_Min_bar;

	vec_XiMin        = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_XiMinBar     = _mat_c_immob_2_xi_immob * conc_bar;

	A = vec_XiMin.segment(_n_xi_Mob + _n_xi_Sorp_bar_li, _n_xi_Min);
	B = vec_XiMinBar.tail(_n_xi_Sorp_bar_ld);

	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde, _n_xi_Min_tilde) = -_vec_XiMinTilde + A - conc_Min_bar - (_mat_Ald * B);
}

// Eq. 3.60
void LocalProblem::residual_xi_Kin(ogsChem::LocalVector & conc_Mob,
								   ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector conc_tmp  = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	conc_tmp        = _mat_c_mob_2_xi_mob * conc_Mob;
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde , _n_xi_Kin)   = - _vec_Xikin + conc_tmp.segment(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin);
}

// Eq. 3.61
void LocalProblem::residual_conc_Sorp(ogsChem::LocalVector & ln_conc_Mob,
									  ogsChem::LocalVector & ln_conc_NonMin_bar,
									  ogsChem::LocalVector & vec_residual)
{

	ogsChem::LocalVector   ln_conc_tmp;
	ln_conc_tmp	 = ogsChem::LocalVector::Zero(_I_mob + _I_NMin_bar);


	//TODO initialize log k values
	ln_conc_tmp.head(_I_mob)      = ln_conc_Mob ;
	ln_conc_tmp.tail(_I_NMin_bar) = ln_conc_NonMin_bar;

	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +_n_xi_Kin, _n_xi_Sorp)  = - _logk_sorp + _mat_Ssorp.transpose() * ln_conc_tmp;
}

// Eq. 3.62
void LocalProblem::residual_conc_Min(ogsChem::LocalVector & ln_conc_Mob,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & conc_Min_bar,
	     	 	 	 	 	 	 	 ogsChem::LocalVector & vec_residual)
{
	size_t        i, idx;
	double        phi(0.0);
	ogsChem::LocalMatrix   mat_S1minT  = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(), _mat_S1min.rows());

	idx  = _n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_xi_Sorp;

	mat_S1minT = _mat_S1min.transpose();

	for (i=0; i < _n_xi_Min; i++)
	{
		phi  = -_logk_min(i) + mat_S1minT.row(i) * ln_conc_Mob;
		vec_residual(idx + i) 	= std::min(phi, conc_Min_bar(i));

	}
}

// Eq. 3.63
void LocalProblem::residual_Eta_bar(ogsChem::LocalVector & conc_NonMin_bar,
									ogsChem::LocalVector & conc_Min_bar,
									ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector      conc_bar  = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);

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

	ogsChem::LocalVector   conc_tmp  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	ogsChem::LocalVector   conc_bar  = ogsChem::LocalVector::Zero(_I_NMin_bar + _n_xi_Min);

	conc_bar.head(_I_NMin_bar)    = conc_NonMin_bar;
	conc_bar.tail(_n_xi_Min)      = conc_Min_bar;

	conc_tmp     = _mat_c_immob_2_xi_immob * conc_bar;
	//XiBarKin is the total mass constrain

	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _n_xi_Kin_bar)  =
			- Xi_Kin_bar + conc_tmp.segment(_n_xi_Sorp_bar + _n_xi_Min_bar, _n_xi_Kin_bar);
	//vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _n_xi_Kin_bar)  = - Xi_Kin_bar + conc_tmp.segment(_n_xi_Sorp_bar + _n_xi_Min_bar, _n_xi_Kin_bar);
}

// Eq. 3.65
void LocalProblem::residual_xi_KinBar_Kin(double deltaT,
										  ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_NonMin_bar,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & Xi_Kin_bar,
										  ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector  conc 		   = ogsChem::LocalVector::Zero(_n_Comp);
	ogsChem::LocalVector  vec_rateKin  = ogsChem::LocalVector::Zero(_J_tot_kin);
	ogsChem::LocalVector  temp_vec1    = ogsChem::LocalVector::Zero(_n_xi_Kin_bar);
	ogsChem::LocalVector  temp_vec2    = ogsChem::LocalVector::Zero(_n_xi_Kin_bar);
	ogsChem::LocalVector  temp_vec3    = ogsChem::LocalVector::Zero(_n_xi_Kin_bar);

	//TODO get theta_water_content
	double theta_water_content (0.5);  //monod2d example

	conc.head	(_I_mob)  					 = conc_Mob;
	conc.segment(_I_mob, _I_NMin_bar) 		 = conc_NonMin_bar;
	conc.tail	(_n_xi_Min_bar) 			 = conc_Min_bar;

	this->reaction_rates(conc, vec_rateKin);

	temp_vec1 = ((theta_water_content * Xi_Kin_bar - theta_water_content * _vec_XiBarKin_old) / deltaT);
	temp_vec2 = (theta_water_content * _mat_A2kin * vec_rateKin);
	temp_vec3 =  temp_vec1 - temp_vec2;
	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _n_xi_Kin_bar) = temp_vec3;

}

//problem specific reaction rates
void LocalProblem::reaction_rates(ogsChem::LocalVector & conc,
								  ogsChem::LocalVector & vec_rateKin)
{
	std::size_t i;
	// then calculate the rates and fill them in the rate vector
	for ( i=0; i < this->_J_tot_kin; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( conc );
		vec_rateKin(i) = this->_list_kin_reactions[i]->getRate();
	}
}

