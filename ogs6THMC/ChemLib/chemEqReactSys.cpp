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
    : _list_eq_reactions(list_eq_reactions), _I(0), _J(0)
{
	// by default, the class is not yet initialized
	isInitialized = false; 

	if ( map_chemComp.size() > 0 && _list_eq_reactions.size() > 0 )
	{   // if there are reactions. 
		// cout how many mobile and how many immobile components
		countComp(map_chemComp); 
        // count how many which type of reactions we have
        countReactions( map_chemComp, _list_eq_reactions ); 
		// make stoichiometric matrix
		buildStoi(map_chemComp, list_eq_reactions);
        // read the logK values
        read_logK(list_eq_reactions); 
        // allocate memory for local residual vector
        _vec_res    = ogsChem::LocalVector::Zero( _I ); 
        // allocate memory for local Jacobi matrix
        _mat_Jacobi = ogsChem::LocalMatrix::Zero( _I, _I ); 
		// flip the initialization flag
		isInitialized = true; 
	}
}

chemEqReactSys::~chemEqReactSys()
{}

void chemEqReactSys::calc_tot_mass(LocalVector & vec_conc_basis, 
                                   LocalVector & vec_conc_second, 
                                   LocalVector & vec_tot_mass)
{
    vec_tot_mass = vec_conc_basis + _matStoi.transpose() * vec_conc_second; 
}

void chemEqReactSys::calc_residual(LocalVector & vec_ln_conc, 
                                   LocalVector & vec_tot_mass_constrain)
{
    size_t i; 
    ogsChem::LocalVector ln_c_basis, ln_c_sec_mob, ln_c_sec_sorp, ln_c_sec_min; 
    ogsChem::LocalVector vec_cur_mass_balance;
    ogsChem::LocalVector vec_conc_basis;
    ogsChem::LocalVector vec_conc_second;
    vec_cur_mass_balance = ogsChem::LocalVector::Zero( _I_basis  ); 
    vec_conc_basis       = ogsChem::LocalVector::Zero( _I_basis  ); 
    vec_conc_second      = ogsChem::LocalVector::Zero( _I_second ); 
    // clean the residual vector
    _vec_res.setZero();
    double res_tmp; 
    size_t res_idx; 
    // residual vector is composed of 3 parts
    // 1) mass action expression of mobile and sorption reactions
    // 2) mass action expression of mineral reactions
    // 3) mass balance expression of basis species
    
    // first split the ln(conc) vector
    ln_c_basis    = vec_ln_conc.head(_I_basis); 
    ln_c_sec_mob  = vec_ln_conc.segment(_I_basis, _I_sec_mob); 
    ln_c_sec_sorp = vec_ln_conc.segment(_I_basis+_I_sec_mob, _I_sec_sorp);  
    ln_c_sec_min  = vec_ln_conc.segment(_I_basis+_I_sec_mob+_I_sec_sorp, _I_sec_min); 
    
    // part 1), _J_mob + _J_sorp  mass action reactions
    for ( i=1; i <  _J_mob ; i++ )
    {
        res_tmp  = -1.0 * _vec_lnK(i) + _matStoi.row(i).dot( ln_c_basis ) - ln_c_sec_mob(i);
        res_idx  = i; 
        _vec_res(res_idx) = res_tmp; 
    }
    for ( i=1; i <  _J_sorp ; i++ )
    {
        res_idx  = i + _J_mob; 
        res_tmp  = -1.0 * _vec_lnK(res_idx) + _matStoi.row(res_idx).dot( ln_c_basis ) - ln_c_sec_sorp(i);
        _vec_res(res_idx) = res_tmp; 
    }

    // part 2), n_react_min mass action reactions
    for ( i=1; i < _J_min; i++ )
    {
        res_idx  = i + _J_mob + _J_sorp; 
        // attention, this is spectial for mineral reactions
        res_tmp  = -1.0 * _vec_lnK(res_idx) + _matStoi.row(res_idx).dot( ln_c_basis ) - ln_c_sec_min(i);
        _vec_res(res_idx) = std::min(res_tmp, exp(ln_c_sec_min(i))); 
    }
    
    // part 3), n_basis mass balance equations
    for ( i=1; i < _I_basis   ; i++ )
        vec_conc_basis(i)                         = exp( ln_c_basis(i)   ); 
    for ( i=1; i < _I_sec_mob ; i++ )
        vec_conc_second(i)                        = exp( ln_c_sec_mob(i) ); 
    for ( i=1; i < _I_sec_sorp; i++ )
        vec_conc_second(_I_sec_mob+i)             = exp( ln_c_sec_sorp(i)); 
    for ( i=1; i < _I_sec_min ; i++ )
        vec_conc_second(_I_sec_mob+_I_sec_sorp+i) = exp( ln_c_sec_min(i) );

    this->calc_tot_mass(vec_conc_basis, 
                        vec_conc_second, 
                        vec_cur_mass_balance); 
    _vec_res.tail(_I_basis) = vec_tot_mass_constrain - vec_cur_mass_balance; 
}

void chemEqReactSys::calc_Jacobi_ana(LocalVector & vec_ln_conc,
                         LocalVector & vec_tot_mass_constrain,
                         LocalVector & vec_res_base)
{
    size_t i, j; 
    const double delta = 1.0e-6; 
    ogsChem::LocalVector ln_c_basis, ln_c_sec_mob, ln_c_sec_sorp, ln_c_sec_min;
    ogsChem::LocalVector vec_unknown_tmp;
    ogsChem::LocalVector vec_conc_tmp; 
    ogsChem::LocalVector mass_balance_res_base; 
    ogsChem::LocalVector mass_balance_res_tmp;
    ogsChem::LocalVector vec_conc_basis;
    ogsChem::LocalVector vec_conc_second;
    vec_conc_tmp          = ogsChem::LocalVector::Zero(_I); 
    mass_balance_res_base = ogsChem::LocalVector::Zero(_I_basis); 
    mass_balance_res_tmp  = ogsChem::LocalVector::Zero(_I_basis); 
    // clean the residual vector
    _mat_Jacobi.setZero(); 

    // first split the ln(conc) vector
    ln_c_basis    = vec_ln_conc.head(_I_basis); 
    ln_c_sec_mob  = vec_ln_conc.segment(_I_basis, _I_sec_mob); 
    ln_c_sec_sorp = vec_ln_conc.segment(_I_basis+_I_sec_mob, _I_sec_sorp);  
    ln_c_sec_min  = vec_ln_conc.segment(_I_basis+_I_sec_mob+_I_sec_sorp, _I_sec_min); 
    
    // fill in the analytical part
    _mat_Jacobi.topLeftCorner (_J, _I_basis) = _matStoi; 
    _mat_Jacobi.topRightCorner(_J, _J).setIdentity(); 
    _mat_Jacobi.topRightCorner(_J, _J) *= -1.0; 
    
    // now the mass balance part
    mass_balance_res_base = vec_res_base.tail(_I_basis);
    // using the numerical increment method to construct Jacobi matrix 
    for ( i=1; i < _I ; i++ )
    {
        // increment the vec_unknown
        vec_unknown_tmp = vec_ln_conc; 
        if ( abs(vec_unknown_tmp(i)) < std::numeric_limits<double>::epsilon() )
        {
            _mat_Jacobi.col(i).tail(_I_basis).setZero(); 
        }
        else
        {
            vec_unknown_tmp(i) += delta * vec_unknown_tmp(i);
            
            for ( j=1; j < _I ; j++ )
                vec_conc_tmp(j) = exp( vec_unknown_tmp(j) ); 
                        
            vec_conc_basis  = vec_conc_tmp.head(_I_basis);
            vec_conc_second = vec_conc_tmp.tail(_I_second);
            calc_tot_mass( vec_conc_basis, vec_conc_second, mass_balance_res_tmp); 
            mass_balance_res_tmp = vec_tot_mass_constrain - mass_balance_res_tmp; 
            _mat_Jacobi.col(i).tail(_I_basis) = ( mass_balance_res_tmp - mass_balance_res_base ) / (delta * vec_unknown_tmp(i));
        }  // end of if else
    }  // end of for i

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Jacobi Matrix: " << std::endl; 
	std::cout << _mat_Jacobi << std::endl;
	// end of debugging-------------------
#endif
}

void chemEqReactSys::buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
	                             std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
{
	size_t i,j, tmp_idx; 
	double tmp_stoi; 
	std::string tmp_str;
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;
	
	// creat the memory for Stoi matrix
	_matStoi_input = LocalMatrix::Zero(_I, _J); 
	
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
				_matStoi_input(tmp_idx,j) = tmp_stoi; 
			}
			else  // this is a basis component
			{
			    _matStoi_input(tmp_idx,j) = 1.0; 
			}

		}  // end of for i
	}  // end of for j

    // number of secondary components
    _I_second = _J; 
    // number of basis components
    _I_basis = _I - _I_second; 
    // number of secondary mobile components
    _I_sec_mob = _I - _I_basis; 
    // organize which components are basis
    // and which are secondary
    _matStoi = _matStoi_input.topRows( _I_basis ).transpose(); 

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Stoichiometric Matrix S: " << std::endl; 
	std::cout << _matStoi << std::endl;
	// end of debugging-------------------
#endif
}

void chemEqReactSys::read_logK(std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions)
{
    size_t i; 
    _vec_lnK = ogsChem::LocalVector::Zero(_J); 
    for ( i=0 ; i < list_eq_reactions.size(); i++ )
        _vec_lnK(i) = list_eq_reactions[i]->get_ln_K(); 
}

void chemEqReactSys::countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
{
	_I_mob = 0; 
    _I_sec_mob = 0; 
	_I_sec_sorp= 0; 
	_I_sec_min = 0; 

	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
	for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
	{
		switch ( it->second->getMobility() )
		{
		case ogsChem::MOBILE: 
			_I_mob++; 
			break;
		case ogsChem::SORPTION: 
			_I_sec_sorp++;
			break;
		case ogsChem::MINERAL: 
			_I_sec_min++;
			break;
		default:
			_I_sec_min++;
			break; 
		}
	}

    _I = _I_mob + _I_sec_sorp + _I_sec_min;
}

void chemEqReactSys::countReactions(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
{
	_J_mob = 0; 
	_J_sorp= 0; 
	_J_min = 0; 

	size_t it; 
	for( it=0; it < list_eq_reactions.size(); it++ )
	{
        switch ( list_eq_reactions[it]->get_type() )
		{
        case ogsChem::MOB_EQ_REACT: 
			_J_mob++; 
			break;
        case ogsChem::SORP_EQ_REACT: 
			_J_sorp++;
			break;
        case ogsChem::MIN_EQ_REACT: 
			_J_min++;
			break;
		default:
			_J_min++;
			break; 
		}  // end of switch
	}  // end of for

    _J = _J_mob + _J_sorp + _J_min;
}



}  // end of namespace