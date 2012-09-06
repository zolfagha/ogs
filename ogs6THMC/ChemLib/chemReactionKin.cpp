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
}

chemReactionKin::~chemReactionKin(void)
{
}

void chemReactionKin::readReactionStr(std::string & reaction_str)
{
	// TODO
	std::string tmp_str;
}

void chemReactionKin::calcReactionRate(void)
{
    // TODO
}

void chemReactionKin::readReactionKRC(ogs5::CKinReact* KRC_reaction)
{
	// local temperary variables
	// size_t tmp_size, i; 

	// go through the list of components and convert it
	// to the local data structure
	

    // copy the stoichiometric vector 
	_vecStoi = KRC_reaction->stochmet;



}

}  // end of namespace