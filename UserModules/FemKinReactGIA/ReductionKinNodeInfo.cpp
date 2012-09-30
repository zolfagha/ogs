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
 
#include "ReductionKinNodeInfo.h"
 
ReductionKinNodeInfo::ReductionKinNodeInfo(size_t node_id, 
	                                       size_t n_comp, 
	                                       size_t n_eta_mob, 
										   size_t n_eta_immob, 
										   size_t n_xi, 
										   ogsChem::chemReductionKin* ReductionKin)
    : _node_id(node_id), _n_comp(n_comp), _n_eta_mob(n_eta_mob), 
	  _n_eta_immob(n_eta_immob), _n_xi(n_xi), _ReductionKin(ReductionKin)
{
	Comp_Conc = DiscreteLib::LocalVector::Zero( _n_comp ); 
	eta_mob   = DiscreteLib::LocalVector::Zero( _n_eta_mob ); ; 
	eta_immob = DiscreteLib::LocalVector::Zero( _n_eta_immob ); ; 
	xi        = DiscreteLib::LocalVector::Zero( _n_xi ); 
}

void ReductionKinNodeInfo::set_comp_conc( size_t comp_idx, double val )
{
	Comp_Conc(comp_idx) = val; 
}

void ReductionKinNodeInfo::transform(void)
{
	_ReductionKin->Conc2EtaXi( Comp_Conc, eta_mob, eta_immob, xi ); 
}

