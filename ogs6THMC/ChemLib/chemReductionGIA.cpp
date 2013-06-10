/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionGIA.h
 *
 * Created on 2013-05-02 by Reza Zolfaghari and Haibing Shao
 */
#include "chemReductionGIA.h"
#include "logog.hpp"
//#define _DEBUG
namespace ogsChem
{
chemReductionGIA::chemReductionGIA(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp,
		                           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions,
		                           std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions)
: _list_eq_reactions(list_eq_reactions), _list_kin_reactions(list_kin_reactions)
{
	// by default, the class is not yet initialized
	isInitialized = false;

	if ( list_chemComp.size() > 0 && (_list_eq_reactions.size() > 0 || _list_kin_reactions.size() > 0) )
	{   // if there are reactions.
		// cout how many mobile and how many immobile components
		countComp(list_chemComp);
		// make stoichiometric matrix
		buildStoi(list_chemComp, list_eq_reactions, list_kin_reactions);
		// calculate the reduction parameters
		update_reductionScheme();
		// get the log k values
		read_logK(list_eq_reactions);
		// flip the initialization flag
		isInitialized = true;
	}
}



void chemReductionGIA::buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp,
								 std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions,
								 std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions)
{
	size_t i,j, tmp_idx;
	double tmp_stoi;
	std::string tmp_str;
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;

    // zero the number of different equilibrium reactions
    _Jmob  = 0; 
    _Jsorp = 0; 
    _Jmin  = 0; 

	// obtain the size info
	_I_tot		 = map_chemComp.size();
	_J_tot_eq 	 = list_eq_reactions.size();
	_J_tot_kin   = list_kin_reactions.size();
	_J_tot		 = _J_tot_eq + _J_tot_kin;

	// creat the memory for Stoi matrix
	_matStoi = LocalMatrix::Zero(_I_tot, _J_tot);

	// based on the reactions, fill in the stoi matrix
	// loop over all reactions
	//Equilibrium reaction
	for ( j=0; j < list_eq_reactions.size() ; j++ )
	{	// for each reaction
        // check which type of equilibrium reaction it is
        if (      list_eq_reactions[j]->get_type() == ogsChem::EqReactType::MOB_EQ_REACT  )
            this->_Jmob++; 
        else if ( list_eq_reactions[j]->get_type() == ogsChem::EqReactType::SORP_EQ_REACT )
            this->_Jsorp++; 
        else if ( list_eq_reactions[j]->get_type() == ogsChem::EqReactType::MIN_EQ_REACT  )
            this->_Jmin++; 

		// find each participating components
		for ( i=0; i <  list_eq_reactions[j]->get_vecCompNames().size() ; i++ ){
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

	//Kinetic reactions
	for ( j=0; j < list_kin_reactions.size() ; j++ )
	{	// for each reaction
		// find each participating components
		for ( i=0; i <  list_kin_reactions[j]->get_vecCompNames().size() ; i++ ){
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

void chemReductionGIA::update_reductionScheme(void)
{
	size_t i,j;

	// divide mobile and immobile parts
	// S1 =S(1:I,:);
	_mat_S1 = _matStoi.topRows( _I_mob );
	// S2 =S(I+1:length(C),:);
	_mat_S2 = _matStoi.bottomRows( _I_NMin_bar + _I_min );

	//preserve the original matrix
	_mat_S1_preserve = _mat_S1;
	_mat_S2_preserve = _mat_S2;

	// prepare the matrix for later extraction
	//extract subspaces of S 1 and 2 matrix
	_mat_S1mob    =_mat_S1.block(0,0,_I_mob,_Jmob);
	_mat_S1sorp   =_mat_S1.block(0,_Jmob,_I_mob,_Jsorp);
	_mat_S1min    =_mat_S1.block(0,_Jmob+_Jsorp,_I_mob,_Jmin);
    _mat_S1kin    =_mat_S1.block(0,_J_tot_eq,_I_mob,_J_tot_kin);

    _mat_S2mob    =_mat_S2.block(0,0,_I_bar,_Jmob);
    _mat_S2sorp   =_mat_S2.block(0,_Jmob,_I_bar,_Jsorp);
	_mat_S2min    =_mat_S2.block(0,_Jmob +_Jsorp,_I_bar,_Jmin);
    _mat_S2kin    =_mat_S2.block(0,_J_tot_eq,_I_bar,_J_tot_kin);

    _Jsorp_li = _Jsorp;
    _mat_S1sorp_li = _mat_S1sorp;
    _matrix_Ssorp  = _matStoi.block(0,_Jmob,_I_mob +_I_NMin_bar,_Jsorp);

    if(_Jsorp)
    {
    	_mat_Ssorp   =_matStoi.block(0,_Jmob,_I_tot,_Jsorp);
    	_mat_Ssorp_li = _mat_Ssorp;   //if _Jsorp = 0

    	// reveal the rank of S1sorp = S1sorpli and S1sorpld
    	Eigen::FullPivLU<LocalMatrix> lu_decomp_Ssorp(_mat_Ssorp);
    	_mat_Ssorp_li = lu_decomp_Ssorp.image(_mat_Ssorp);

    	//_mat_Ssorp_li = _mat_Ssorp_ast;
    	//Linear_Indep =[];
    	//Linear_Dep   =[];
    	for (j= 0; j < _Jsorp; j++)
    	{
    		for (i=0; i< _I_tot; i++)
    		{
    			if (_mat_Ssorp_li(i,j) != _mat_Ssorp(i,j))
    			{
    				_mat_Ssorp_ld (i,j) = _mat_Ssorp(i,j);
    			}
    		}
    	}

    	_mat_S1sorp_li.resize(0,0);   //clean the allocated memory
    	_mat_S1sorp_li = _mat_Ssorp_li.topRows(_I_mob);
    	_mat_S2sorp_li = _mat_Ssorp_li.bottomRows(_I_bar);
    	_Jsorp_li   = _mat_Ssorp_li.cols();
    	_Jsorp_ld   = _mat_Ssorp_ld.cols();

    	_mat_S1sorp = _mat_S1sorp_li; // if _Jsorp_ld = 0
    	_mat_S2sorp = _mat_S2sorp_li; // if _Jsorp_ld = 0

    	// construct the global _mat_Ssorp matrix, right parts contains the linearly independent and left part contains the linearly dependent reactions.
    	_matrix_Ssorp.leftCols(_Jsorp_li)  = _mat_Ssorp_li;
    	_matrix_Ssorp.rightCols(_Jsorp_ld) = _mat_Ssorp_ld;


    	if( _mat_Ssorp_ld.cols())
    	{
    		_mat_S1sorp_ld = _mat_Ssorp_ld.topRows(_I_mob);
    		_mat_S2sorp_ld = _mat_Ssorp_ld.bottomRows(_I_bar);

    		_mat_S1sorp << _mat_S1sorp_li,_mat_S1sorp_ld;
    		_mat_S2sorp << _mat_S2sorp_li,_mat_S2sorp_ld;

    	}
    }

    //reconstruct the original matrix S1 based on the new arrangement (li,ld)
	_mat_S1.block(0,0,_I_mob,_Jmob) = _mat_S1mob;
	_mat_S1.block(0,_Jmob,_I_mob,_Jsorp) = _mat_S1sorp;
	_mat_S1.block(0,_Jmob+_Jsorp,_I_mob,_Jmin) = _mat_S1min;
    _mat_S1.block(0,_J_tot_eq,_I_mob,_J_tot_kin) = _mat_S1kin;

    _mat_S2.block(0,0,_I_bar,_Jmob) = _mat_S2mob;
    _mat_S2.block(0,_Jmob,_I_bar,_Jsorp) = _mat_S2sorp;
	_mat_S2.block(0,_Jmob+_Jsorp,_I_bar,_Jmin) = _mat_S2min;
    _mat_S2.block(0,_J_tot_eq,_I_bar,_J_tot_kin) = _mat_S2kin;

    _mat_S1kin_ast = _mat_S1kin;
    _mat_S2kin_ast = _mat_S2kin;
    _Jkin_ast      = this->_J_tot_kin;  // HS: is this correct?

    if( _J_tot_kin > 0 )
    {
    	//clean the memory
    	_mat_S1kin_ast.resize(0,0);
    	_mat_S2kin_ast.resize(0,0);

    	// reveal the rank of S1kin
    	Eigen::FullPivLU<LocalMatrix> lu_decomp_Skin(_mat_Skin);
    	_mat_Skin_li = lu_decomp_Skin.image(_mat_Skin);
        _mat_S1kin_ast = _mat_Skin_li.topRows(_I_mob);
        _mat_S2kin_ast = _mat_Skin_li.bottomRows(_I_bar);
        _Jkin_ast      = _mat_Skin_li.cols();
    }

	// creat the memory for Stoi matrix
	_mat_S1_ast = LocalMatrix::Zero(_I_mob, _Jmob + _Jsorp_li + _Jmin +_Jkin_ast);
	_mat_S2_ast = LocalMatrix::Zero(_I_bar, _Jsorp + _Jmin + _Jkin_ast);

	// construct S1* and S2* which consists of max linearly independent columns of S1 and 2
    _mat_S1_ast.block(0,0,_I_mob,_Jmob) 							= _mat_S1mob;
    _mat_S1_ast.block(0,_Jmob,_I_mob,_Jsorp_li) 					= _mat_S1sorp_li;
    _mat_S1_ast.block(0,_Jmob + _Jsorp_li,_I_mob,_Jmin) 			= _mat_S1min;
    _mat_S1_ast.block(0,_Jmob + _Jsorp_li + _Jmin,_I_mob,_Jkin_ast) = _mat_S1kin_ast;

    _mat_S2_ast.block(0,0,_I_bar,_Jsorp) 						  = _mat_S2sorp;
    _mat_S2_ast.block(0,_Jsorp, _I_bar ,_Jmin) 		   	  = _mat_S2min;
    _mat_S2_ast.block(0,_Jsorp + _Jmin,_I_bar,_Jkin_ast) = _mat_S2kin_ast;

    //the location of these matrices should not be changed. Cut the zero parts
	_mat_S2sorp = _mat_S2sorp.topRows(_I_NMin_bar);
	_mat_S2kin_ast = _mat_S2kin_ast.topRows(_I_NMin_bar);

#ifdef _DEBUG
    std::cout << "_mat_S1: "    << std::endl;
	std::cout << _mat_S1 << std::endl;
	std::cout << "_mat_S2: "    << std::endl;
	std::cout << _mat_S2 << std::endl;
	std::cout << "_mat_S1_ast: "    << std::endl;
	std::cout << _mat_S1_ast << std::endl;
	std::cout << "_mat_S2_ast: "    << std::endl;
	std::cout << _mat_S2_ast << std::endl;
#endif

	// Calculate the s1T = S_i^T matrix consisting of a max set of linearly
	// independent columns that are orthogonal to each column of S_i^*.
	// s1T = orthcomp(s1);
	_mat_S1_orth = orthcomp( _mat_S1_ast );
    // s2T = orthcomp(s2);
	_mat_S2_orth = orthcomp( _mat_S2_ast );

#ifdef _DEBUG
	std::cout << "_mat_S1_orth: "    << std::endl;
	std::cout << _mat_S1_orth << std::endl;
	std::cout << "_mat_S2_orth: "    << std::endl;
	std::cout << _mat_S2_orth << std::endl;
#endif

	_mat_c_mob_2_eta_mob     = ( _mat_S1_orth.transpose() * _mat_S1_orth ).fullPivHouseholderQr().solve(_mat_S1_orth.transpose());
    _mat_c_immob_2_eta_immob = ( _mat_S2_orth.transpose() * _mat_S2_orth ).fullPivHouseholderQr().solve(_mat_S2_orth.transpose());
	_mat_c_mob_2_xi_mob      = ( _mat_S1_ast.transpose()  * _mat_S1_ast  ).fullPivHouseholderQr().solve(_mat_S1_ast.transpose());
	_mat_c_immob_2_xi_immob  = ( _mat_S2_ast.transpose()  * _mat_S2_ast  ).fullPivHouseholderQr().solve(_mat_S2_ast.transpose());

#ifdef _DEBUG
	std::cout << "_mat_c_mob_2_eta_mob: "    << std::endl;
	std::cout << _mat_c_mob_2_eta_mob << std::endl;
	std::cout << "_mat_c_immob_2_eta_immob: "    << std::endl;
	std::cout << _mat_c_immob_2_eta_immob << std::endl;
	std::cout << "_mat_c_mob_2_xi_mob: "    << std::endl;
	std::cout << _mat_c_mob_2_xi_mob << std::endl;
	std::cout << "_mat_c_immob_2_xi_immob: "    << std::endl;
	std::cout << _mat_c_immob_2_xi_immob << std::endl;
#endif
/*
 * known indeces : Jmob, Jsorp, Jmin, Jkin, IbarNmin
 */


	// now store the corresponding size
	this->_n_eta             = _mat_S1_orth.cols();
	this->_n_eta_bar         = _mat_S2_orth.cols();
	this->_n_xi_mobile       = _mat_S1_ast.cols();
	this->_n_xi_immobile     = _mat_S2_ast.cols();

	//_li and _ast is max num of linearly independent columns
    _Jsorp_li = _mat_Ssorp_li.cols();
    _Jsorp_ld = _mat_Ssorp_ld.cols();
    _Jkin_ast = _mat_S1kin_ast.cols();
    _J_tot_li = _Jmob + _Jsorp_li + _Jmin + _Jkin_ast;

	//known
	this->_n_xi_Mob          = _Jmob;
	this->_n_xi_Sorp_bar     = _Jsorp;
	this->_n_xi_Sorp         = _Jsorp_li;
	this->_n_xi_Min		     = _Jmin;
	this->_n_xi_Min_tilde    = _Jmin;
	this->_n_xi_Min_bar      = _Jmin;

	//calculated
	this->_n_xi_Sorp_tilde   = _Jsorp_li;
	this->_n_xi_Sorp_bar_ld  = _Jsorp_ld;
	this->_n_xi_Sorp_bar_li  = _Jsorp_li;
	this->_n_xi_Kin          = _Jkin_ast;
	this->_n_xi_Kin_bar      = _Jkin_ast;

	this->_n_xi_global       = _Jsorp_li + _Jsorp_li + _Jmin + _Jmin + _Jkin_ast;
	this->_n_xi_local        = _Jmob + _Jsorp + _Jmin + _Jkin_ast;
	this->_Jeq_li            = _Jmob + _Jsorp_li + _Jmin;


	// now calculated A1 and A2
	_mat_A1 = ( _mat_S1_ast.transpose() * _mat_S1_ast ).fullPivHouseholderQr().solve(_mat_S1_ast.transpose()) * _mat_S1 ;
	_mat_A2 = ( _mat_S2_ast.transpose() * _mat_S2_ast ).fullPivHouseholderQr().solve(_mat_S2_ast.transpose()) * _mat_S2 ;

#ifdef _DEBUG
	std::cout << "A_1: "    << std::endl;
	std::cout << _mat_A1 << std::endl;
	std::cout << "A_2: "    << std::endl;
	std::cout << _mat_A2 << std::endl;
#endif

	//extract the A subspaces
	_mat_A1sorp       = _mat_A1.block(_Jmob,_J_tot_eq,_Jsorp_li,_J_tot_kin);
	_mat_A1min        = _mat_A1.block(_Jmob+_Jsorp_li,_J_tot_eq,_Jmin,_J_tot_kin);
	_mat_A1kin        = _mat_A1.block(_Jmob+_Jsorp_li+_Jmin,_J_tot_eq,_Jkin_ast,_J_tot_kin);
	_mat_Ald	      = _mat_A1.block(_Jmob+_Jsorp_li,_Jmob+_Jsorp_li,_Jmin,_Jsorp_ld);
	_mat_A2sorp       = _mat_A2.block(0,_Jmob+_Jsorp+_Jmin,_Jsorp,_J_tot_kin);
	_mat_A2sorpli     = _mat_A2sorp.block(0,0,_Jsorp_li,_J_tot_kin);
	_mat_A2sorpld     = _mat_A2sorp.block(_Jsorp_li,0,_Jsorp_ld,_J_tot_kin);
	_mat_A2kin        = _mat_A2.block(_Jsorp+_Jmin,_Jmob+_Jsorp+_Jmin,_Jkin_ast,_J_tot_kin);

}


LocalMatrix chemReductionGIA::orthcomp( LocalMatrix & inMat )
{
	// initialize it so that they have the same dimension
    LocalMatrix outMat;

	Eigen::FullPivLU<LocalMatrix> lu_decomp(inMat.transpose());

	outMat = lu_decomp.kernel();

	// notice that if the kernel returns a matrix with dimension zero,
	// then the returned matrix will be a column-vector filled with zeros
	// therefore we do a safety check here, and set the number of columns to zero
	if ( outMat.cols() == 1 && outMat.col(0).norm() == 0.0 )
		outMat = LocalMatrix::Zero(inMat.rows(), 0);

	return outMat;
}


void chemReductionGIA::countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
{
	_I_mob = 0;
	_I_NMin_bar= 0;
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
			_I_NMin_bar++;
			break;
		case ogsChem::MINERAL:
			_I_min++;
			break;
		default:
			_I_min++;
			break;
		}
	}
    _I_bar = _I_NMin_bar + _I_min;
}

void chemReductionGIA::read_logK(std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions)
{
    size_t i;
    ogsChem::LocalVector vec_logK;
    vec_logK = ogsChem::LocalVector::Zero(_Jmob + _Jsorp + _Jmin);

    for ( i=0 ; i < list_eq_reactions.size(); i++ )
        vec_logK(i) = list_eq_reactions[i]->get_ln_K();

    _logk_mob  = vec_logK.head(_Jmob);
    _logk_sorp = vec_logK.segment(_Jmob, _Jsorp);
    _logk_min  = vec_logK.tail(_Jmin);
}

void chemReductionGIA::Conc2EtaXi(ogsChem::LocalVector &local_conc,
								  ogsChem::LocalVector &local_eta,
								  ogsChem::LocalVector &local_eta_bar,
								  ogsChem::LocalVector &local_xi_global,
								  ogsChem::LocalVector &local_xi_local)
{
	// declare local temp variable
	ogsChem::LocalVector local_c_mob, local_c_immob;
	ogsChem::LocalVector local_xi, local_xi_bar;
	ogsChem::LocalVector local_xi_Sorp_bar_li, local_xi_Sorp_bar_ld;

	// divide c1 and c2
	local_c_mob   = local_conc.topRows(    this->_I_mob );
	local_c_immob = local_conc.bottomRows( this->_I_NMin_bar + this->_I_min );

	// convert eta_mob and xi_mob
	local_eta   = _mat_c_mob_2_eta_mob * local_c_mob;
	local_xi    = _mat_c_mob_2_xi_mob  * local_c_mob;

	// convert eta_immob and xi_immob
	local_eta_bar = _mat_c_immob_2_eta_immob * local_c_immob;
    local_xi_bar  = _mat_c_immob_2_xi_immob  * local_c_immob;
    //0 means it is a vector, one column.
    local_xi_Mob  = local_xi.segment( 0,this->_n_xi_Mob);
    local_xi_Sorp = local_xi.segment( this->_n_xi_Mob,this->_n_xi_Sorp);
    local_xi_Min  = local_xi.segment( this->_n_xi_Mob + this->_n_xi_Sorp,this->_n_xi_Min);
    local_xi_Kin  =	local_xi.segment( this->_n_xi_Mob + this->_n_xi_Sorp + this->_n_xi_Min,this->_n_xi_Kin);

    local_xi_Sorp_bar =  local_xi_bar.segment( 0,this->_n_xi_Sorp);
    local_xi_Min_bar  =  local_xi_bar.segment( this->_n_xi_Sorp_bar,this->_n_xi_Min_bar);
    local_xi_Kin_bar  =  local_xi_bar.segment( this->_n_xi_Sorp_bar + this->_n_xi_Min_bar,this->_n_xi_Kin_bar);


    local_xi_Sorp_bar_li = local_xi_Sorp_bar.topRows(this->_Jsorp_li);
    local_xi_Sorp_bar_ld = local_xi_Sorp_bar.bottomRows(this->_Jsorp_ld);
    local_xi_Sorp_tilde  = local_xi_Sorp - local_xi_Sorp_bar.topRows(this->_Jsorp_li);
    local_xi_Min_tilde   = local_xi_Min  - local_xi_Min_bar - (_mat_Ald * local_xi_Sorp_bar_ld);

    local_xi_local.segment( 0,this->_n_xi_Mob) 				     = local_xi_Mob;
    local_xi_local.segment(this->_n_xi_Mob,this->_n_xi_Sorp_bar) = local_xi_Sorp_bar;
    local_xi_local.segment(this->_n_xi_Mob + this->_n_xi_Sorp_bar,this->_n_xi_Min_bar) = local_xi_Min_bar;
    local_xi_local.segment(this->_n_xi_Mob + this->_n_xi_Sorp_bar + this->_n_xi_Min_bar,this->_n_xi_Kin_bar) = local_xi_Kin_bar;

    local_xi_global.segment( 0,this->_n_xi_Sorp_tilde) 														     = local_xi_Sorp_tilde;
    local_xi_global.segment( this->_n_xi_Sorp_tilde,this->_n_xi_Min_tilde) 				  					     = local_xi_Min_tilde;
    local_xi_global.segment( this->_n_xi_Sorp_tilde + this->_n_xi_Min_tilde,this->_n_xi_Sorp) 				     = local_xi_Sorp;
    local_xi_global.segment( this->_n_xi_Sorp_tilde + this->_n_xi_Min_tilde + this->_n_xi_Sorp,this->_n_xi_Min) 				   = local_xi_Min;
    local_xi_global.segment( this->_n_xi_Sorp_tilde + this->_n_xi_Min_tilde + this->_n_xi_Sorp + this->_n_xi_Min,this->_n_xi_Kin)  = local_xi_Kin;

}

void chemReductionGIA::EtaXi2Conc(ogsChem::LocalVector &local_eta,
	                              ogsChem::LocalVector &local_eta_bar,
								  ogsChem::LocalVector &local_xi_global,
								  ogsChem::LocalVector &local_xi_local,
								  ogsChem::LocalVector &local_conc )
{
	// declare local temp variable
	ogsChem::LocalVector local_xi, local_xi_bar, local_c_mob, local_c_immob;
	local_xi 	  = ogsChem::LocalVector::Zero(_n_xi_Mob +_n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	local_c_mob   = ogsChem::LocalVector::Zero(_I_mob);
	local_xi_bar  = ogsChem::LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	local_c_immob = ogsChem::LocalVector::Zero(_I_NMin_bar + _I_min);

	//local_xi_Mob  = local_xi_local.segment( 0,this->_n_xi_Mob);
	//local_xi_Sorp = local_xi_global.segment(this->_n_xi_Sorp_tilde + this->_n_xi_Min,this->_n_xi_Sorp);
	//local_xi_Min  = local_xi_global.segment(this->_n_xi_Sorp_tilde + this->_n_xi_Min + this->_n_xi_Sorp,this->_n_xi_Min);
	//local_xi_Kin  = local_xi_global.segment(this->_n_xi_Sorp_tilde + this->_n_xi_Min + this->_n_xi_Sorp+ this->_n_xi_Min,this->_n_xi_Kin);

	//local_xi_Sorp_bar = local_xi_local.segment(this->_n_xi_Mob,this->_n_xi_Sorp_bar);
	//local_xi_Min_bar  = local_xi_local.segment(this->_n_xi_Mob+this->_n_xi_Sorp_bar,this->_n_xi_Min_bar);
	//local_xi_Kin_bar  = local_xi_local.segment(this->_n_xi_Mob+this->_n_xi_Sorp_bar+this->_n_xi_Min_bar,this->_n_xi_Kin_bar);

	local_xi.head	(this->_n_xi_Mob) 					  				   = local_xi_local.segment( 0,this->_n_xi_Mob);;
	local_xi.segment(this->_n_xi_Mob, this->_n_xi_Sorp) 		           = local_xi_global.segment(this->_n_xi_Sorp + this->_n_xi_Min,this->_n_xi_Sorp);;
	local_xi.segment(this->_n_xi_Mob + this->_n_xi_Sorp, this->_n_xi_Min ) = local_xi_global.segment(this->_n_xi_Sorp + this->_n_xi_Min + this->_n_xi_Sorp,this->_n_xi_Min);
	local_xi.tail	(this->_n_xi_Kin) 									   = local_xi_global.segment(this->_n_xi_Sorp + this->_n_xi_Min + this->_n_xi_Sorp+ this->_n_xi_Min,this->_n_xi_Kin);

	local_xi_bar.head	(_n_xi_Sorp_bar) 					 =  local_xi_local.segment(this->_n_xi_Mob,this->_n_xi_Sorp_bar);
	local_xi_bar.segment(_n_xi_Sorp_bar, _n_xi_Min_bar) 	 =  local_xi_local.segment(this->_n_xi_Mob+this->_n_xi_Sorp_bar,this->_n_xi_Min_bar);
	local_xi_bar.tail	(_n_xi_Kin_bar) 					 =  local_xi_local.segment(this->_n_xi_Mob+this->_n_xi_Sorp_bar+this->_n_xi_Min_bar,this->_n_xi_Kin_bar);

	local_c_mob   = _mat_S1_ast * local_xi   + _mat_S1_orth * local_eta;
	local_c_immob = _mat_S2_ast * local_xi_bar + _mat_S2_orth * local_eta_bar;

	local_conc.topRows( this->_I_mob ) = local_c_mob;
	local_conc.bottomRows( this->_I_NMin_bar + this->_I_min ) = local_c_immob;

    // testing if the non-negative stablilization will help?
    for (size_t i=0; i < local_conc.rows(); i++)    {
        if ( local_conc(i) < 0.0 )
            local_conc(i) = 1.0e-99;
    }
    // end of testing

}



void chemReductionGIA::Calc_Kin_Rate(ogsChem::LocalVector &local_eta,
	                                ogsChem::LocalVector &local_eta_bar,
									ogsChem::LocalVector &local_xi_global,
									ogsChem::LocalVector &local_xi_local,
									ogsChem::LocalVector &local_rate_vec)
{
	size_t i;

	// the size of vec_rates is equal to the number of kinetic equations
	ogsChem::LocalVector vec_rates = ogsChem::LocalVector::Zero(_J_tot_kin);
	// the local temp concentration vector
	ogsChem::LocalVector vec_conc = ogsChem::LocalVector::Zero(_I_tot);

	// first convert these eta and xi to concentrations
	EtaXi2Conc( local_eta,
		        local_eta_bar,
				local_xi_global,
				local_xi_local,
				vec_conc );

	// then calculate the rates and fill them in the rate vector
	for ( i=0; i < this->_J_tot_kin; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( vec_conc );
		vec_rates(i) = this->_list_kin_reactions[i]->getRate();
	}

	//Note: since A sub spaces are different for different equations, it should be multiplied later
	// multiply the rate vector with the A matrix to get rate for xi_mob and xi_immob
	//xi_mob_rate   = _matA1 * vec_rates;
	//xi_immob_rate = _matA2 * vec_rates;
}


}  // end of namespace

