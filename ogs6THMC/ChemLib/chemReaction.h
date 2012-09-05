/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReaction.h
 *
 * Created on 2012-08-23 by Haibing Shao
 */
#ifndef CHEMREACTION_H
#define CHEMREACTION_H

#include "chemcomp.h"

namespace ogsChem
{

class chemReaction
{
public:
	/**
      * constructor of the class
      */
	chemReaction(void);
	/**
      * destructor of the class
      */
    ~chemReaction(void); 

protected:
	/**
      * virtual class, read from a reaction string 
      */
	virtual void readReactionStr(std::string & reaction_str); 

	/**
      * set the reactions 
      * @param pcomp is the pointer to an instance of chemical component
      * @param stoi  is the stoichiometric coefficient of coresponding component
      */
    void addComp( ChemComp* pComp, double stoi ); 
    
	/**
      * vector of components in this reaction
      */
    std::vector<ChemComp*> _vecComponents; 

    /**
      * vector of stoichiometric coefficient before the component
      */
    std::vector<double> _vecStoi; 

};

} // end of namespace


#endif