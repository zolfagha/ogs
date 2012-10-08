/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Local_ODE_Xi_immob.h
 *
 * Created on 2012-10-08 by Haibing Shao
 */
 
#ifndef LOCAL_ODE_XI_IMMOB_H
#define LOCAL_ODE_XI_IMMOB_H

#include "ChemLib/chemReductionKin.h"

class Local_ODE_Xi_immob
{
public: 
	typedef MathLib::LocalVector MyVector; 
    Local_ODE_Xi_immob(ogsChem::chemReductionKin* reductionKin)
		: _reductionKin(reductionKin)
    {
		_vec_eta_mob      = MyVector::Zero( _reductionKin->get_n_eta_mob() ); 
	    _vec_eta_immob    = MyVector::Zero( _reductionKin->get_n_eta_immob() ); 
	    _vec_xi_mob       = MyVector::Zero( _reductionKin->get_n_xi_mob() ); 
	    _vec_xi_immob     = MyVector::Zero( _reductionKin->get_n_xi_immob() ); 
		_vec_dxi_immob_dt = MyVector::Zero( _reductionKin->get_n_xi_immob() ); 	 
	}
	 
	~Local_ODE_Xi_immob()
    {
        _reductionKin = NULL; 
    }

	void update_eta_xi( MyVector & vec_eta_mob, 
		                MyVector & vec_eta_immob, 
						MyVector & vec_xi_mob, 
						MyVector & vec_xi_immob  )
	{
		_vec_eta_mob   = vec_eta_mob; 
		_vec_eta_immob = vec_eta_immob; 
	    _vec_xi_mob    = vec_xi_mob; 
	    _vec_xi_immob  = vec_xi_immob;
	}

	MyVector operator() (double & time, MyVector & y )
	{
	    _vec_xi_immob = y; 
		this->_reductionKin->Calc_Xi_immob_Rate( _vec_eta_mob, _vec_eta_immob, _vec_xi_mob, _vec_xi_immob, _vec_dxi_immob_dt); 
		return _vec_dxi_immob_dt; 
	}


private:     
	ogsChem::chemReductionKin* _reductionKin; 
	MyVector _vec_eta_mob; 
	MyVector _vec_eta_immob; 
	MyVector _vec_xi_mob; 
	MyVector _vec_xi_immob;
	MyVector _vec_dxi_immob_dt; 
 
}; 
 
 #endif