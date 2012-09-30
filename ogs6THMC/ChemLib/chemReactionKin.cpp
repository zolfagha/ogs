/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReaction.h
 *
 * Created on 2012-09-05 by Haibing Shao
 */
 
#include "chemReactionKin.h"
 
namespace ogsChem
{

chemReactionKin::chemReactionKin(void)
{
	// set initial values
	_rate = 0.0;
	
	_rate_constant = 0.0; 
    
	_rate_constant_order - 0.0; 
}

chemReactionKin::~chemReactionKin(void)
{
}

void chemReactionKin::readReactionStr(std::string & reaction_str)
{
	// TODO
	std::string tmp_str;
}

void chemReactionKin::calcReactionRate(ogsChem::LocalVector & vec_Comp_Conc)
{
	if ( this->_kinReactType == KinReactType::Monod )
		this->_rate = calcReactionRateMonod(vec_Comp_Conc); 
	
}

double chemReactionKin::calcReactionRateMonod(ogsChem::LocalVector & vec_Comp_Conc)
{
	size_t i, comp_idx; 
	double rate = 1.0; 
	double conc; 

	// loop over all the monod term
	for (i=0; i<this->_vec_Monod_Comps_Idx.size(); i++ )
	{
		comp_idx = _vec_Monod_Comps_Idx[i]; 
		conc = vec_Comp_Conc( comp_idx ); 
		// rate *= [ c / ( c + k_c) ]^order
		rate *= pow( conc / (conc + this->_vec_Monod_Comps_Conc[i]), this->_vec_Monod_Comps_order[i] ); 
	}  // end of for i
		
	// retrun rate
	return rate; 
}

void chemReactionKin::readReactionKRC(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, ogs5::CKinReact* KRC_reaction)
{
	// get the list of components 
	_vecCompNames = KRC_reaction->reactionpartner; 

	// find the components based on component names
	size_t i; 
	std::string str_comp; 
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator comp_iterator; 
	for (i=0; i<_vecCompNames.size(); i++)
	{
		str_comp = _vecCompNames[i];
		comp_iterator = list_chemComp.find(str_comp);

		if ( comp_iterator != list_chemComp.end() )
		{
			ogsChem::ChemComp* component = comp_iterator->second; 
			this->_vecComponents.push_back(component);		
		}  // end of if	
	}  // end of for i

    // copy the stoichiometric vector 
	_vecStoi = KRC_reaction->stochmet;

	// read the rate parameters
	if ( KRC_reaction->getType() == "monod" )
	{
		this->_kinReactType = KinReactType::Monod; 

		// loop over the monod term, 
		for (size_t i=0; i < KRC_reaction->monod.size(); i++ )
		{
			size_t monodComp_idx; 
			double monod_comp_conc(0.0), monod_term_order(1.0); 
			// read the corresponding components
			std::string comp_str = KRC_reaction->monod[i]->species; 
			// find this component
			for (size_t j=0; j < this->_vecComponents.size(); j++ )
			{
				if ( this->_vecComponents[j]->get_name() == comp_str )
				{
					monodComp_idx = this->_vecComponents[j]->getIndex(); 
				    break; 
				}
			}  // end of for j

			// read the monod component concentration
			monod_comp_conc = KRC_reaction->monod[i]->concentration; 
			// read the monod term order
			monod_term_order = KRC_reaction->monod[i]->order; 

			this->_vec_Monod_Comps_Idx.push_back(monodComp_idx);
			this->_vec_Monod_Comps_Conc.push_back( monod_comp_conc ); 
			this->_vec_Monod_Comps_order.push_back( monod_term_order ); 
		}  // end of for i
	}  // end of if KRC_reaction
	else
	{
		this->_kinReactType = KinReactType::NoType; 
	}

	// read the rate constant
	this->_rate_constant = KRC_reaction->rateconstant; 
	// read the rate order
	this->_rate_constant_order = KRC_reaction->rateorder;


}




}  // end of namespace