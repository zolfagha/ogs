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
						 size_t n_xi, 
						 ogsChem::chemReductionKin* ReductionKin); 
	~ReductionKinNodeInfo(); 

	void set_comp_conc( size_t comp_idx, double val); 

	size_t get_node_id(void) { return _node_id; }

	void transform(void); 

protected: 


private: 
	const size_t _node_id; 
	const size_t _n_comp; 
	const size_t _n_eta_mob; 
	const size_t _n_eta_immob; 
	const size_t _n_xi; 

	DiscreteLib::LocalVector Comp_Conc; 
	DiscreteLib::LocalVector eta_mob; 
	DiscreteLib::LocalVector eta_immob; 
	DiscreteLib::LocalVector xi; 

	ogsChem::chemReductionKin* _ReductionKin;
}; 

#endif  // end of ifndef