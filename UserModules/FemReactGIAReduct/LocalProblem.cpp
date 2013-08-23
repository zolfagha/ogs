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


void LocalProblem::solve_LocalProblem_Newton_LineSearch(ogsChem::LocalVector & x,
														ogsChem::LocalVector & vec_tot_mass_constrain,
														std::size_t & node_idx ,
														double deltaT,
														const double iter_tol,
														const double rel_tol,
														const double max_iter)
{
    ogsChem::LocalVector x_new, vec_residual, vec_AI;
    ogsChem::LocalVector dx;
    ogsChem::LocalVector ln_conc_Mob, ln_conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar;
//TODO initialize water content and deltaT
    // number of iterations
    size_t j, iter, n_unknowns;
    const double alpha (0.5);
    double d_norm, d1_norm;

    n_unknowns     = _n_Comp + _n_xi_Kin_bar;
    //x_new          = ogsChem::LocalVector::Ones( n_unknowns );
    x_new          = ogsChem::LocalVector::Zero( n_unknowns );
    dx             = ogsChem::LocalVector::Ones( n_unknowns );
    //dx             = ogsChem::LocalVector::Zero( n_unknowns );
    _mat_Jacobian  = ogsChem::LocalMatrix::Zero(n_unknowns, n_unknowns);
    //vec_residual   = ogsChem::LocalVector::Ones( n_unknowns );
    vec_residual   = ogsChem::LocalVector::Zero( n_unknowns );
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
    ln_conc_Mob     	= x.head(_I_mob );
    //NonMin refers to non mineral (sorbed + kinetic) concentrations.
    ln_conc_NonMin_bar  = x.segment(_I_mob, _I_NMin_bar );
    conc_Min_bar     	= x.segment(_I_mob + _I_NMin_bar, _I_min );
    // xi kin bar is both unknown and mass constrain!
    Xi_Kin_bar       	= x.tail(_n_xi_Kin_bar);


    // start solving the system
    iter = 0; 

    // now updating the saturation index and minerals
    this->update_minerals_conc_AI( x, vec_AI );

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "x Vector: \n";
	//std::cout << x << std::endl;

	//std::cout << "vec_AI Vector: \n";
	//std::cout << vec_AI << std::endl;
	// end of debugging-------------------
#endif


    // evaluate the residual
    this->calc_residual(x, vec_residual);

    // evaluate norm of residual vector
    d_norm = vec_residual.norm();

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "Residual Vector: \n";
	//std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

	// save the previous values
	x_new  =  x;

    while (true)
    {
        #ifdef _DEBUG
            // display the residual
//             std::cout << "Iteration #" << iter << "||res|| = " << d_norm << "||delta_x|| = " << dx.norm() << std::endl;
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
        this->calc_Jacobian(x, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "_mat_Jacobian: \n";
	//std::cout << _mat_Jacobian << std::endl;
	// end of debugging-------------------
#endif

        // solving for increment
        this->Solv_Minimization( node_idx, _mat_Jacobian, vec_residual, dx );

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "dx Vector: \n";
	//std::cout << dx << std::endl;
	// end of debugging-------------------
#endif


        // increment of unknowns
        this->increment_unknown( x, dx, x_new ); 

#ifdef _DEBUG
	// debugging--------------------------
	// std::cout << "x_new Vector: \n";
	// std::cout << x_new << std::endl;
	// end of debugging-------------------
#endif

    	// updating the saturation index and minerals
		this->update_minerals_conc_AI( x_new, vec_AI );

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "x_new: \n";
	//std::cout << x_new << std::endl;
	// end of debugging-------------------
#endif

        // evaluate residual with x_new
        this->calc_residual(x_new, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
	//std::cout << "x_new Vector: \n";
	//std::cout << x_new << std::endl;

	//std::cout << "vec_residual: \n";
	//std::cout << vec_residual << std::endl;
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
            this->update_minerals_conc_AI( x_new, vec_AI );
            // evaluate residual with x_new
            this->calc_residual(x_new, vec_residual);

		#ifdef _DEBUG
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

void LocalProblem::calc_residual(ogsChem::LocalVector & vec_unknowns,
								 ogsChem::LocalVector & vec_residual)
{
	ogsChem::LocalVector ln_conc_Mob, ln_conc_NonMin_bar, conc_Mob, conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar;
	conc_Mob 		    = ogsChem::LocalVector::Zero(_I_mob);
	conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);

    ln_conc_Mob         = vec_unknowns.segment(0, _I_mob );
    //NonMin refers to sorbed and kinetic concentrations.
    ln_conc_NonMin_bar  = vec_unknowns.segment(_I_mob, _I_NMin_bar );
    conc_Min_bar        = vec_unknowns.segment(_I_mob + _I_NMin_bar, _I_min );
    Xi_Kin_bar   	    = vec_unknowns.segment(_I_mob + _I_NMin_bar + _I_min, _n_xi_Kin_bar );

    // convert the ln mobile conc to mobile conc
     this->cal_exp_conc_vec(_I_mob, ln_conc_Mob, conc_Mob);
     // convert the ln immobile conc to immobile conc
     this->cal_exp_conc_vec(_I_NMin_bar, ln_conc_NonMin_bar, conc_Min_bar);

    // Eq. 3.55
    this->residual_conc_Mob			(ln_conc_Mob, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "residual_conc_Mob: \n";
//	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.56
    this->residual_Eta				(conc_Mob, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "residual_Eta: \n";
//	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.57 - 58
    this->residual_xi_Sorp_tilde	(conc_Mob, conc_NonMin_bar, conc_Min_bar, vec_residual);
    // Eq. 3.59
    this->residual_xi_Min_tilde		(conc_Mob, conc_NonMin_bar, conc_Min_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
//	std::cout << "residual_xi_Min_tilde: \n";
    //	std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.60
    this->residual_xi_Kin			(conc_Mob, vec_residual);
    // Eq. 3.61
    this->residual_conc_Sorp		(ln_conc_Mob, ln_conc_NonMin_bar, vec_residual);
    // Eq. 3.62
    this->residual_conc_Min			(ln_conc_Mob, conc_Min_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
    //	std::cout << "residual_conc_Min: \n";
    //std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.63
    this->residual_Eta_bar			(conc_NonMin_bar, conc_Min_bar, vec_residual);

#ifdef _DEBUG
	// debugging--------------------------
    //std::cout << "residual_Eta_bar: \n";
    //std::cout << vec_residual << std::endl;
	// end of debugging-------------------
#endif

    // Eq. 3.64
    this->residual_xi_KinBar_Eq		(conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar, vec_residual);
    // Eq. 3.65
    this->residual_xi_KinBar_Kin	(conc_Mob, conc_NonMin_bar, conc_Min_bar, Xi_Kin_bar, vec_residual);

}  // end of function calc_residual

void LocalProblem::calc_Jacobian(ogsChem::LocalVector & vec_x,
//							     ogsChem::LocalVector & vec_AI,
							     ogsChem::LocalVector & vec_residual)
{
	const double delta_xi = 1.0e-8;
    std::size_t i;
    ogsChem::LocalVector vec_x_incremented, vec_residual_incremented;
    vec_residual_incremented = vec_residual;

	for (i = 0; i < _mat_Jacobian.cols(); i++)
	{
		vec_x_incremented  = vec_x;

//		// numerical protection
//		if( vec_x_incremented.norm() < 1.0e-16)
		{
			vec_x_incremented(i) += delta_xi * vec_x_incremented.norm();
			this->calc_residual(vec_x_incremented, vec_residual_incremented);
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
    ogsChem::LocalVector ln_conc_Mob, ln_conc_NonMin_bar, conc_Min_bar;
    ogsChem::LocalMatrix mat_S1min_transposed;
    //logConc_Mob = ogsChem::LocalVector::Zero(_I_mob);
    mat_S1min_transposed = ogsChem::LocalMatrix::Zero(_mat_S1min.cols(),_mat_S1min.rows());
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
    ogsChem::LocalVector temp, conc_bar, Conc_Mob, conc_NonMin_bar;
    ogsChem::LocalMatrix xi_sorp_bar_Ald;
    ogsChem::LocalVector xi_mobile;
    temp 				= ogsChem::LocalVector::Zero(_n_xi_Kin_bar + _n_xi_Sorp_bar + _n_xi_Min_bar);
    conc_bar 			= ogsChem::LocalVector::Zero(_n_xi_Kin_bar + _n_xi_Sorp_bar + _n_xi_Min_bar);
    Conc_Mob 			= ogsChem::LocalVector::Zero(_I_mob);
    conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
    xi_mobile           = ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
    xi_sorp_bar_Ald     = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

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
	ogsChem::LocalVector                 vec_XiSorp;
	ogsChem::LocalVector                 vec_XiSorpBarLI;
	ogsChem::LocalVector                 conc_bar;
	conc_bar = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);

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
	ogsChem::LocalVector                 vec_XiMin;
	ogsChem::LocalVector                 vec_XiMinBar;
	ogsChem::LocalVector                 conc_bar, A, B;
	conc_bar 		= ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);
	vec_XiMin 		= ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	vec_XiMinBar 	= ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	A 				= ogsChem::LocalVector::Zero(_n_xi_Min);
	B 				= ogsChem::LocalVector::Zero(_n_xi_Sorp_bar_ld);

	conc_bar.head(_n_xi_Sorp_bar)    = conc_NonMin_bar;
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
//	ogsChem::LocalVector   conc_tmp;
//	conc_tmp 		= ogsChem::LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
//
//	conc_tmp        = _mat_c_mob_2_xi_mob * conc_Mob;
//	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde , _n_xi_Kin)   = - _vec_Xikin + conc_tmp.segment(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin);
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
	double        phi;
	ogsChem::LocalMatrix   mat_S1minT;
	mat_S1minT 	= ogsChem::LocalMatrix::Zero(_mat_S1min.cols(), _mat_S1min.rows());

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

//	ogsChem::LocalVector   conc_tmp;
//	ogsChem::LocalVector   conc_bar;
//	conc_bar  = ogsChem::LocalVector::Zero(_I_NMin_bar + _n_xi_Min);
//	conc_tmp  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
//
//	conc_bar.head(_I_NMin_bar)    = conc_NonMin_bar;
//	conc_bar.tail(_n_xi_Min)         = conc_Min_bar;
//
//	conc_tmp     = _mat_c_immob_2_xi_immob * conc_bar;
//	//XiBarKin is the total mass constrain
//	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar, _n_xi_Kin_bar)  = - Xi_Kin_bar + conc_tmp.segment(_n_xi_Sorp_bar + _n_xi_Min_bar, _n_xi_Kin_bar);
}

// Eq. 3.65
void LocalProblem::residual_xi_KinBar_Kin(ogsChem::LocalVector & conc_Mob,
										  ogsChem::LocalVector & conc_NonMin_bar,
										  ogsChem::LocalVector & conc_Min_bar,
										  ogsChem::LocalVector & Xi_Kin_bar,
										  ogsChem::LocalVector & vec_residual)
{
//	ogsChem::LocalVector   conc, vec_rateKin;
//	conc 		 = ogsChem::LocalVector::Zero(_n_Comp);
//	vec_rateKin  = ogsChem::LocalVector::Zero(_n_xi_Kin);
//
//	//TODO get theta_waterContent
//	double theta_waterContent (0.3);
//
//	conc.head	(_I_mob)  					 = conc_Mob;
//	conc.segment(_I_mob, _I_NMin_bar) 		 = conc_NonMin_bar;
//	conc.tail	(_n_xi_Min_bar) 			 = conc_Min_bar;
//
//	this->reaction_rates(conc, vec_rateKin);
//
//	vec_residual.segment(_n_xi_Mob + _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde +  _n_xi_Kin + _n_xi_Sorp + _n_xi_Min + _n_eta_bar + _n_xi_Kin_bar, _n_xi_Kin_bar)
//			      = ((theta_waterContent * Xi_Kin_bar - (theta_waterContent * _vec_XiBarKin_old)) / deltaT) - (theta_waterContent * _mat_A2kin * vec_rateKin);
}

//problem specific reaction rates
void LocalProblem::reaction_rates(ogsChem::LocalVector & conc,
								  ogsChem::LocalVector & vec_rateKin)
{
//double sub_A, sub_B, sub_C, biomass, umax, k_subA, k_subB, kdec;
//ogsChem::LocalVector R;
//ogsChem::LocalVector::Zero();
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
