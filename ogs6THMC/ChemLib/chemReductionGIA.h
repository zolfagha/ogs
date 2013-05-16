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
#ifndef CHEMREDUCTIONGIA_H
#define CHEMREDUCTIONGIA_H

#include "chemconst.h"
#include "chemReactionKin.h"
#include "chemReactionEq.h"
#include "chemcomp.h"
#include "BaseLib/OrderedMap.h"

namespace ogsChem
{

class chemReductionGIA
{
public:
	/**
      * constructor of the class
      */
	chemReductionGIA(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp,
		             std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions,
		             std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions);
	
	/**
      * destructor of the class
      */
	~chemReductionGIA(void);

    /**
      * convert concentration vector to eta and xi vector

	void Conc2EtaXi(ogsChem::LocalVector &local_conc,
					ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalVector &local_xi_Mob,
					ogsChem::LocalVector &local_xi_Sorp_tilde,
					ogsChem::LocalVector &local_xi_Sorp_bar,
					ogsChem::LocalVector &local_xi_Min_tilde,
					ogsChem::LocalVector &local_xi_Min_bar,
					ogsChem::LocalVector &local_xi_Kin,
					ogsChem::LocalVector &local_xi_Kin_bar);
*/
	void Conc2EtaXi(ogsChem::LocalVector &local_conc,
					ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
				    ogsChem::LocalMatrix &local_xi_global,
					ogsChem::LocalMatrix &local_xi_local);

	/**
      * convert eta and xi vectors to concentration vector

	void EtaXi2Conc(ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalVector &local_xi_Mob,
					ogsChem::LocalVector &local_xi_Sorp_tilde,
					ogsChem::LocalVector &local_xi_Sorp_bar,
					ogsChem::LocalVector &local_xi_Min_tilde,
					ogsChem::LocalVector &local_xi_Min_bar,
					ogsChem::LocalVector &local_xi_Kin,
					ogsChem::LocalVector &local_xi_Kin_bar,
					ogsChem::LocalVector &local_conc);
 */
	void EtaXi2Conc(ogsChem::LocalVector &local_eta,
					ogsChem::LocalVector &local_eta_bar,
					ogsChem::LocalMatrix &local_xi_global,
					ogsChem::LocalMatrix &local_xi_local,
					ogsChem::LocalVector &local_conc);
    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 

	/**
	 * calculate the reaction rates
	 */
	void Calc_Kin_Rate(ogsChem::LocalVector &local_eta,
		               ogsChem::LocalVector &local_eta_bar,
					   ogsChem::LocalMatrix &local_xi_global,
					   ogsChem::LocalMatrix &local_xi_local,
					   ogsChem::LocalVector &local_rate_vec);

	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I_tot; };
	
	/**
      * get the length of eta_mob
      */
	size_t get_n_eta_mob(void) {return _n_eta;};

	/**
      * get the length of eta_immob
      */
	size_t get_n_eta_immob(void) {return _n_eta_bar;};

	/**
      * get the length of eta_mob
      */
	size_t get_n_xi_mobile(void) {return _n_xi_mobile;};

	/**
      * get the length of eta_mob
      */
	size_t get_n_xi_immobile(void) {return _n_xi_immobile;};

	/**
      * get the length of xi_Mob = Jmob
      */
	size_t get_n_xi_Mob(void) {return _n_xi_Mob; };

	/**
      * get the length of xi_Sorp_tilde = Jsorpli
      */
	size_t get_n_xi_Sorp(void) {return _n_xi_Sorp; };

	/**
      * get the length of xi_Min_bar = Jmin
      */
	size_t get_n_xi_Min(void) {return _n_xi_Min; };


	/**
      * get the length of xi_Sorp_tilde = Jsorpli
      */
	size_t get_n_xi_Sorp_tilde(void) {return _n_xi_Sorp_tilde; };

	/**
      * get the length of xi_Sorp_bar = Jsorp
      */
	size_t get_n_xi_Sorp_bar(void) {return _n_xi_Sorp_bar; };

	/**
      * get the length of xi_Min_tilde = Jmin
      */
	size_t get_n_Min_tilde(void) {return _n_xi_Min_tilde; };

	/**
      * get the length of xi_Min_bar = Jmin
      */
	size_t get_n_xi_Min_bar(void) {return _n_xi_Min_bar; };
	/**
      * get the length of xi_ki = j1kin
      */
	size_t get_n_xi_Kin(void) {return _n_xi_Kin; };

	/**
      * get the length of xi_Kin_bar = j2kin
      */
	size_t get_n_xi_Kin_bar(void) {return _n_xi_Kin_bar; };

private:
	/**
      * private flag indicating initialization
      */
	bool isInitialized; 

	/**
      * a list of all kinetic / equilibrium reactions
      */
	std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions;

	/**
      * a list of all equilibrium reactions
      */
	std::vector<ogsChem::chemReactionEq*> & _list_eq_reactions;

    /**
      * stoichiometric matrix S
      */	
	LocalMatrix _matStoi; 

	/**
      * sub-matrixes of S, 
      */	
	LocalMatrix _mat_S1, _mat_S1mob,_mat_S2mob, _mat_Ssorp,_mat_Ssorp_ast,_mat_Ssorp_li,_mat_Ssorp_ld, _mat_S1sorp_li, _mat_S1sorp_ld, _mat_S1sorp, _mat_S1min, _mat_S2min, _mat_S2, _mat_S2sorp, _mat_S1kin, _mat_S2kin,
	            _mat_Skin,_mat_Skin_li, _mat_S2sorp_li,_mat_S2sorp_ld,_mat_S1_preserve, _mat_S2_preserve;

	/**
      * LocalVector sub-vectors of xi and eta,
      */
	LocalMatrix local_xi_Mob, local_xi_Sorp,local_xi_Sorp_li,local_xi_Sorp_ld,local_xi_Sorp_tilde,local_xi_Sorp_bar,local_xi_Min,local_xi_Min_tilde,local_xi_Min_bar,local_xi_Kin,local_xi_Kin_bar,
	            local_xi,local_xi_bar;

	LocalVector local_conc,local_eta,local_eta_bar;


	/**
      * Max linearly independent columns (astride/ast) of the matrixes
      */
	LocalMatrix _mat_S1kin_ast, _mat_S2kin_ast, _mat_S1_ast, _mat_S2_ast;

	/**
      * orthogonal complementary of the matrixes
      */
	LocalMatrix _mat_S1_orth, _mat_S2_orth, _mat_S2_orth_tilde;

	/**
      * matrix for the calculation of right hand side rate transformation
      */
	LocalMatrix _mat_A1,_mat_A2,_mat_A1sorp, _mat_A1min, _mat_A1kin, _mat_A2kin, _mat_Ald, _mat_A2sorp, _mat_A2sorpli, _mat_A2sorpld;

	/**
      * transformation matrix for the calculation between eta, xi and concentrations
      */
	LocalMatrix _mat_c_mob_2_eta_mob,
		        _mat_c_immob_2_eta_immob,   
	            _mat_c_mob_2_xi_mob, 
				_mat_c_immob_2_xi_immob;

	/**
      * the size of eta and xi vector
      */
	size_t _n_eta,_n_eta_bar, _n_xi_mobile, _n_xi_immobile, _n_xiLocal,  _n_xiGlobal, _n_xi_Mob, _n_xi_Sorp_tilde,_n_xi_Sorp, _n_xi_Sorp_bar,_n_xi_Min, _n_xi_Min_tilde, _n_xi_Min_bar, _n_xi_Kin, _n_xi_Kin_bar;

	/**
      * _I_tot is the number of all components and _J_tot is the number of all reactions in the system
      */
	//size_t _I_tot, _J_tot;

	/**
      * _I_tot is the number of all components and _J_tot is the number of all reactions in the system
      */
	//size_t 	Linear_Indep, Linear_Dep;
	/**
      * _I_tot is the number of all components and _J_tot is the number of all reactions in the system
      */
	size_t _J_tot_eq, _J_tot_kin,_J_tot_li,_J_tot;

	/**
      * number of mobile, sorption and mineral and kinetic reactions
      */
	size_t _Jeq_li, _Jmob, _Jsorp, _Jsorp_li, _Jsorp_ld, _Jmin, _Jkin, _Jkin_ast;

	/**
      * number of mobile and immobile components and species
      */
	size_t _I_tot,_I_mob, _I_sorp,_I_min, _I_immob, _I_immob_min, _I_immob_nonMin;

	/**
      * construct stoichiometric matrix out of list of components and reactions
      */
	void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & list_chemComp, 
		           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions,
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
