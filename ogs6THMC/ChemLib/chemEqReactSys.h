/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemEqReactSys.h
 *
 * Created on 2013-03-28 by Haibing Shao
 */
#ifndef CHEM_EQ_REACT_SYS_H
#define CHEM_EQ_REACT_SYS_H

#include "chemconst.h"
#include "chemcomp.h"
#include "chemReactionEq.h"
#include "BaseLib/OrderedMap.h"

namespace ogsChem
{

class chemEqReactSys
{
public:
	/**
      * constructor of the class
      */
	chemEqReactSys(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
                   std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions);
	
	/**
      * destructor of the class
      */
	~chemEqReactSys(void); 

    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 
	
	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I; }; 
	
private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

	/**
      * a list of all kinetic reactions
      */
	std::vector<ogsChem::chemReactionEq*> & _list_eq_reactions; 

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

	/**
      * _I is the number of components and _J is the number of reactions
      */
	size_t _I, _J; 

	/**
      * number of mobile, sorption and mineral reactions
      */
	size_t _J_mob, _J_sorp, _J_min;

    /**
      * number of mobile, sorption and mineral components
      */
	size_t _I_mob, _I_sorp, _I_min;

	/**
      * construct stoichiometric matrix out of list of components and reactions
      */
	void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions);

	/**
      * count how many mobile and immobile components
      */
	void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp); 

};

}

#endif