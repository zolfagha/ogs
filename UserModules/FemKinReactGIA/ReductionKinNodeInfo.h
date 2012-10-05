/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepKinReduction.h
 *
 * Created on 2012-09-27 by Haibing Shao
 */

#ifndef REDUCTION_KIN_NODE_INFO_H
#define REDUCTION_KIN_NODE_INFO_H

#include "ChemLib/chemReductionKin.h"

class ReductionKinNodeInfo
{
public: 
	ReductionKinNodeInfo(size_t node_id, 
		                 size_t n_comp, 
	                     size_t n_eta_mob, 
						 size_t n_eta_immob, 
						 size_t n_xi_mob,
						 size_t n_xi_immob,
						 ogsChem::chemReductionKin* ReductionKin); 
	~ReductionKinNodeInfo(); 

	void set_comp_conc( size_t comp_idx, double val); 

	size_t get_node_id(void) { return _node_id; }

	double get_eta_mob_value( size_t eta_mob_idx ) { return _eta_mob(eta_mob_idx); }

	double get_eta_immob_value( size_t eta_immob_idx ) { return _eta_immob(eta_immob_idx); }

	double get_xi_mob_value( size_t xi_mob_idx ) { return _xi_mob(xi_mob_idx); }

	double get_xi_immob_value( size_t xi_immob_idx ) { return _xi_immob(xi_immob_idx); }

	void transform(void); 

protected: 


private: 
	const size_t _node_id; 
	const size_t _n_comp; 
	const size_t _n_eta_mob; 
	const size_t _n_eta_immob; 
	const size_t _n_xi_mob; 
	const size_t _n_xi_immob; 

	MathLib::LocalVector _Comp_Conc; 
	MathLib::LocalVector _eta_mob; 
	MathLib::LocalVector _eta_immob; 
	MathLib::LocalVector _xi_mob; 
	MathLib::LocalVector _xi_immob; 

	ogsChem::chemReductionKin* _ReductionKin;
}; 

#endif  // end of ifndef