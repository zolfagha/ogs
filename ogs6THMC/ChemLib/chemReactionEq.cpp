/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReactionEq.h
 *
 * Created on 2013-03-18 by Haibing Shao
 */
 
#include "chemReactionEq.h"
#include "logog.hpp"

namespace ogsChem
{

chemReactionEq::chemReactionEq()
{
	// set initial values
	_eq_const_K = 1.0000001;
	_ln_K = log(_eq_const_K);  

}

chemReactionEq::~chemReactionEq(void)
{
}

void chemReactionEq::readReactionStr(std::string & reaction_str)
{
	// TODO
	std::string tmp_str;
}

void chemReactionEq::readReactionKRC(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                                      ogs5::CKinReact* KRC_reaction)
{
//	// get the list of components 
//	_vecCompNames = KRC_reaction->reactionpartner; 
//
//	// find the components based on component names
//	size_t i; 
//	std::string str_comp; 
//	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator comp_iterator; 
//	for (i=0; i<_vecCompNames.size(); i++)
//	{
//		str_comp = _vecCompNames[i];
//		comp_iterator = list_chemComp.find(str_comp);
//
//		if ( comp_iterator != list_chemComp.end() )
//		{
//			ogsChem::ChemComp* component = comp_iterator->second; 
//			this->_vecComponents.push_back(component);		
//		}  // end of if	
//	}  // end of for i
//
//    // copy the stoichiometric vector 
//	_vecStoi = KRC_reaction->stochmet;
//
//	// read the rate parameters
//	if ( KRC_reaction->getType() == "monod" )
//	{
//		this->_kinReactType = ogsChem::Monod; 
//
//        // find the bacteria component
//        std::string str_bacteria = KRC_reaction->bacteria_name; 
//        this->_idx_bacteria = list_chemComp.find( str_bacteria )->second->getIndex(); 
//
//        if ( this->_idx_bacteria >= list_chemComp.size() )
//        {
//            ERR("When reading Monod reaction, bacteria name not found! ");
//            exit(1);
//        }
//
//		// loop over the monod term, 
//		for (size_t i=0; i < KRC_reaction->monod.size(); i++ )
//		{
//			size_t monod_comp_idx; 
//			double monod_comp_conc(0.0), monod_term_order(1.0); 
//			// read the corresponding components
//			std::string comp_str = KRC_reaction->monod[i]->species; 
//			// find this component
//            // get its index value in all component list
//            monod_comp_idx = list_chemComp.find( comp_str )->second->getIndex();
//			
//            if ( monod_comp_idx >= list_chemComp.size() )
//            {
//                ERR("When reading Monod terms, component name not found! ");
//                exit(1);
//            }
//
//			// read the monod component concentration
//			monod_comp_conc = KRC_reaction->monod[i]->concentration; 
//			// read the monod term order
//			monod_term_order = KRC_reaction->monod[i]->order; 
//
//			this->_vec_Monod_Comps_Idx.push_back(   monod_comp_idx   );
//			this->_vec_Monod_Comps_Conc.push_back(  monod_comp_conc  ); 
//			this->_vec_Monod_Comps_order.push_back( monod_term_order ); 
//		}  // end of for i
//	}  // end of if KRC_reaction
//	else
//	{
//		this->_kinReactType = ogsChem::NoType; 
//	}
//
//	// read the rate constant
//	this->_rate_constant = KRC_reaction->rateconstant; 
//	// read the rate order
//	this->_rate_constant_order = KRC_reaction->rateorder;
//    // read the decay rate
//    this->_decay_rate = KRC_reaction->decay_rate; 
}


}  // end of namespace