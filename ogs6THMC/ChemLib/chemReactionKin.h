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
 
#ifndef CHEM_REACTION_KIN_H
#define CHEM_REACTION_KIN_H

#include "chemReaction.h"
#include "ogsFileIO\FemIO\ogs5\rf_kinreact.h"

namespace ogsChem
{

class chemReactionKin : public chemReaction
{
public: 
	/**
      * constructor and destructor
      */
	chemReactionKin(void);
	~chemReactionKin(void);

	/**
      * force the program to calculate the reaction rate
      */
	void calcReactionRate(void); 

	/**
      * return the rate of current reaction. 
      */
	double getRate(void) { return _rate; }

	/**
      * override reading str function
      */
	void readReactionStr(std::string & reaction_str); 

	/**
      * reading one reaction from KRC data structure
      */
	void readReactionKRC(ogs5::CKinReact* KRC_reaction); 

private:
	/**
      * the rate of reaction
      */
	double _rate; 
};

}  // end of namespace

#endif