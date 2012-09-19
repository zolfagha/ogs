/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionKin.cpp
 *
 * Created on 2012-08-23 by Haibing Shao
 */

#include "chemReductionKin.h"
#include "logog.hpp"

namespace ogsChem
{

chemReductionKin::chemReductionKin(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
	                               std::vector<ogsChem::chemReactionKin*> & list_kin_reactions)
{
	// by default, the class is not yet initialized
	isInitialized = false; 

	if ( map_chemComp.size() > 0 && list_kin_reactions.size() > 0 )
	{   // if there are reactions. 
		// cout how many mobile and how many immobile components
		countComp(map_chemComp); 
		// make stoichiometric matrix
		buildStoi(map_chemComp, list_kin_reactions); 
		// calculate the reduction parameters	
		update_reductionScheme(); 
		// flip the initialization flag
		isInitialized = true; 
	}
}

chemReductionKin::~chemReductionKin()
{}

void chemReductionKin::buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
	                             std::vector<ogsChem::chemReactionKin*> & list_kin_reactions)
{
	size_t i,j, tmp_idx; 
	double tmp_stoi; 
	std::string tmp_str;
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;
	
	// obtain the size info
	_I = map_chemComp.size(); 
	_J = list_kin_reactions.size(); 
	// creat the memory for Stoi matrix
	_matStoi = LocalMatrix::Zero(_I, _J); 
	
	// based on the reactions, fill in the stoi matrix
	// loop over all reactions
	for ( j=0; j < list_kin_reactions.size(); j++ )
	{	// for each reaction
		// find each participating components 
		for ( i=0; i < list_kin_reactions[j]->get_vecCompNames().size(); i++ ){
			tmp_str  = list_kin_reactions[j]->get_vecCompNames()[i]; 
			tmp_Comp = map_chemComp.find(tmp_str);
			tmp_idx  = tmp_Comp->second->getIndex(); 
			if ( list_kin_reactions[j]->get_vecStoi().size() > 2 )
			{
				// normal reactions
				tmp_stoi = list_kin_reactions[j]->get_vecStoi()[i];
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

void chemReductionKin::update_reductionScheme(void)
{
	// devide mobile and immobile parts
	// S1 =S(1:I,:);
	_matS_1 = _matStoi.topRows( _I_mob );
	// S2 =S(I+1:length(C),:);
	_matS_2 = _matStoi.bottomRows( _I_min ); 

	// reveal the rank of S1 and S2
	Eigen::FullPivLU<LocalMatrix> lu_decomp_S1(_matS_1); 
	_mat_s_1 = lu_decomp_S1.image(_matS_1); 
	Eigen::FullPivLU<LocalMatrix> lu_decomp_S2(_matS_2); 
	_mat_s_2 = lu_decomp_S2.image(_matS_2); 

#ifdef _DEBUG
	std::cout << "_mat_s_1: "    << std::endl; 
	std::cout << _mat_s_1 << std::endl;
	std::cout << "_mat_s_2: "    << std::endl; 
	std::cout << _mat_s_2 << std::endl;
#endif

	// Calculate the s1T = S_i^T matrix consisting of a max set of linearly
	// independent columns that are orthogonal to each column of S_i^*.
	// s1T = orthcomp(s1);
	_matS_1_ast = orthcomp( _mat_s_1 );
    // s2T = orthcomp(s2);
	_matS_2_ast = orthcomp( _mat_s_2 );

#ifdef _DEBUG
	std::cout << "s_1*: "    << std::endl; 
	std::cout << _matS_1_ast << std::endl;
	std::cout << "s_2*: "    << std::endl; 
	std::cout << _matS_2_ast << std::endl;
#endif

	// now store the corresponding size
	this->_n_eta_mob   = _mat_s_1.cols(); 
	this->_n_eta_immob = _mat_s_2.cols();
	this->_n_xi_mob    = _matS_1_ast.cols(); 
	this->_n_xi_immob  = _matS_2_ast.cols();
	this->_n_eta       = _n_eta_mob + _n_eta_immob; 
	this->_n_xi        = _n_xi_mob  + _n_xi_immob;
}

LocalMatrix chemReductionKin::orthcomp( LocalMatrix & inMat )
{
	// initialize it so that they have the same dimension
    LocalMatrix outMat; 

	Eigen::FullPivLU<LocalMatrix> lu_decomp(inMat.transpose());

	outMat = lu_decomp.kernel(); 

	return outMat;
}

void chemReductionKin::countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
{
	_I_mob = 0; 
	_I_sorp= 0; 
	_I_min = 0; 

	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
	for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
	{
		switch ( it->second->getMobility() )
		{
		case ogsChem::Comp_Mobility::MOBILE: 
			_I_mob++; 
			break;
		case ogsChem::Comp_Mobility::SORPTION: 
			_I_sorp++;
			break;
		case ogsChem::Comp_Mobility::MINERAL: 
			_I_min++;
			break;
		default:
			_I_min++;
			break; 
		}
	}
}

}  // end of namespace