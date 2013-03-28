/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemEqReactSys.cpp
 *
 * Created on 2013-03-26 by Haibing Shao
 */

#include "chemEqReactSys.h"
#include "logog.hpp"

namespace ogsChem
{

chemEqReactSys::chemEqReactSys(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
                               std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
    : _list_eq_reactions(list_eq_reactions)
{
	// by default, the class is not yet initialized
	isInitialized = false; 

	if ( map_chemComp.size() > 0 && _list_eq_reactions.size() > 0 )
	{   // if there are reactions. 
		// cout how many mobile and how many immobile components
		countComp(map_chemComp); 
		// make stoichiometric matrix
		buildStoi(map_chemComp, list_eq_reactions); 
		// flip the initialization flag
		isInitialized = true; 
	}
}

chemEqReactSys::~chemEqReactSys()
{}

void chemEqReactSys::buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
	                             std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
{
	size_t i,j, tmp_idx; 
	double tmp_stoi; 
	std::string tmp_str;
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;
	
	// obtain the size info
	_I = map_chemComp.size(); 
	_J = list_eq_reactions.size(); 
	// creat the memory for Stoi matrix
	_matStoi = LocalMatrix::Zero(_I, _J); 
	
	// based on the reactions, fill in the stoi matrix
	// loop over all reactions
	for ( j=0; j < list_eq_reactions.size(); j++ )
	{	// for each reaction
		// find each participating components 
		for ( i=0; i < list_eq_reactions[j]->get_vecCompNames().size(); i++ ){
			tmp_str  = list_eq_reactions[j]->get_vecCompNames()[i]; 
			tmp_Comp = map_chemComp.find(tmp_str);
			tmp_idx  = tmp_Comp->second->getIndex(); 
			if ( list_eq_reactions[j]->get_vecStoi().size() > 2 )
			{
				// normal reactions
				tmp_stoi = list_eq_reactions[j]->get_vecStoi()[i];
				// and put them into Stoi matrix
				_matStoi(tmp_idx,j) = tmp_stoi; 
			}
			else  // this is a basis component
			{
			    _matStoi(tmp_idx,j) = 1.0; 
			}

		}  // end of for i
	}  // end of for j

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Stoichiometric Matrix S: " << std::endl; 
	std::cout << _matStoi << std::endl;
	// end of debugging-------------------
#endif
}


void chemEqReactSys::countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
{
	_I_mob = 0; 
	_I_sorp= 0; 
	_I_min = 0; 

	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
	for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
	{
		switch ( it->second->getMobility() )
		{
		case ogsChem::MOBILE: 
			_I_mob++; 
			break;
		case ogsChem::SORPTION: 
			_I_sorp++;
			break;
		case ogsChem::MINERAL: 
			_I_min++;
			break;
		default:
			_I_min++;
			break; 
		}
	}
}




}  // end of namespace