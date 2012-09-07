/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionKin.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */
#ifndef CHEMREDUCTIONKIN_H
#define CHEMREDUCTIONKIN_H

#include "chemconst.h"
#include "chemReactionKin.h"
#include "chemcomp.h"
#include "BaseLib/OrderedMap.h"

namespace ogsChem
{

class chemReductionKin
{
public:
	
	/**
      * constructor of the class
      */
	chemReductionKin(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		             std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions);
	
	/**
      * destructor of the class
      */
	~chemReductionKin(void); 

    /**
      * convert concentration vector to eta and xi vector
      */
	void Conc2EtaXi();

	/**
      * convert eta and xi vector to concentration vector
      */
	void EtaXi2Conc(); 

    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 

private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

	/**
      * sub-matrixes of S
      */	
	LocalMatrix _matS_1, _matS_2, _matS_1_bar, _matS_2_bar; 

	/**
      * complementary orthorgonal matrixes
      */
	LocalMatrix _matS_1_ast, _matS_2_ast, _matS_1_bar_ast, _matS_2_bar_ast; 

	/**
      * _I is the number of components and _J is the number of reactions
      */
	size_t _I, _J; 

	/**
      * number of mobile, sorption and mineral components
      */
	size_t _I_mob, _I_sorp, _I_min;

	/**
      * construct stoichiometric matrix out of list of components and reactions
      */
	void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		           std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions);

	/**
      * calculate the intemediate parameters for the reduction scheme
      */
	void update_reductionScheme(void); 

	/**
      * return the orthogonal complement of the given matrix
      */
	LocalMatrix orthcomp( LocalMatrix & inMat ); 

	/**
      * count how many mobile and immobile components
      */
	void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp); 

};

}

#endif