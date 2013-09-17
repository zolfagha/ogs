#ifndef CHEM_EQ_REACT_SYS_ACTIVITY_H
#define CHEM_EQ_REACT_SYS_ACTIVITY_H

#include "chemconst.h"
#include "chemcomp.h"
#include "chemReactionEq.h"
#include "BaseLib/OrderedMap.h"
#include "chemActivityModelUnity.h"

namespace ogsChem
{

/**
  * The main difference btw. this class and the class chemEqReactSys is, 
  * we us [logC; c_min_bar] as the unknown vector. 
  * The size of this vector is equal to the number of components. 
  * by keeping the log, we will keep concentration values to positive
  */
class chemEqReactSysActivity
{
public:
    /**
      * constructor of the class
      */
	chemEqReactSysActivity(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
                           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions, 
                           ogsChem::chemActivityModelAbstract                  *a )
      : _list_eq_reactions(list_eq_reactions), _I(0), _J(0), _activity_model(a)
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
            _vec_res    = ogsChem::LocalVector::Zero( _I_basis + _I_sec_min ); 
            // allocate memory for local Jacobi matrix
            _mat_Jacobi = ogsChem::LocalMatrix::Zero( _I, _I ); 
		    // allocate memory for vector AI
            if ( this->_I_sec_min > 0 )
                _AI = ogsChem::LocalVector::Zero(_I_sec_min); 
            // flip the initialization flag
            isInitialized = true; 
	    }
    };
    
    /**
      * destructor of the class
      */
    ~chemEqReactSys(void){}; 
    /**
      * whether the reduction scheme has been initialized
      */
	bool IsInitialized(void) {return isInitialized;}; 
	
	/**
      * get the number of components
      */
	size_t get_n_Comp(void) {return _I; }; 

    /**
      * get the number of mobile components
      */
	size_t get_n_Comp_mob(void) {return _I_mob; }; 

    /**
      * calculate the total mass using 
      * the concentration vector of basis, and
      * the concentration vector of secondary variables
      */
    void calc_tot_mass(LocalVector & vec_conc_basis, 
                       LocalVector & vec_conc_second, 
                       LocalVector & vec_tot_mass)
    {
        vec_tot_mass = vec_conc_basis - _matStoi.transpose() * vec_conc_second; 
    }; 

    /**
      * calculate the residual of the reaction system using 
      * the concentration vector of all components and
      * the total concentration constrain of the basis species
      */
    void calc_residual(LocalVector & vec_unknowns, 
                       LocalVector & vec_tot_mass_constrain, 
                       LocalVector & vec_residual)
    {
        /*
        size_t i; 
        double res_tmp, phi;
        ogsChem::LocalVector c_basis, c_sec_min, c_second; 
        ogsChem::LocalVector ln_c_basis, ln_c_sec_mob, ln_c_sec_sorp; 
        ogsChem::LocalVector vec_conc_basis, vec_conc; 
        ogsChem::LocalVector vec_cur_mass_balance;
        ogsChem::LocalVector lnK_min;
        ogsChem::LocalMatrix Stoi_mob, Stoi_sorp, Stoi_min; 

        vec_conc             = ogsChem::LocalVector::Zero( _I ); 
        vec_cur_mass_balance = ogsChem::LocalVector::Zero( _I_basis  ); 
        ln_c_basis           = ogsChem::LocalVector::Zero( _I_basis );
        vec_conc_basis       = ogsChem::LocalVector::Zero( _I_basis  ); 
        c_second             = ogsChem::LocalVector::Zero( _I_second );

        // now updating the saturation index and minerals
        // this->update_AI( vec_unknowns ); 
        // this->update_minerals( vec_unknowns, vec_tot_mass_constrain );

        // now split the unknown vector
        c_basis    = vec_unknowns.head(_I_basis); 
        for (i=0; i < (size_t)c_basis.rows(); i++)
            ln_c_basis(i) = std::log(c_basis(i));
        c_sec_min  = vec_unknowns.tail(_I_sec_min); 

        // part 0), calculate the concentration of secondary 
        // non-mineral components
        Stoi_mob  = _matStoi.topRows(    _J_mob );
        Stoi_sorp = _matStoi.middleRows( _J_mob, _J_sorp ); 
        Stoi_min  = _matStoi.bottomRows( _J_min );
    
        // TODO, do the activity correction

        // calculate the secondary mobile component concentrations
        ln_c_sec_mob  = _vec_lnK.head( _J_mob ) - Stoi_mob * ln_c_basis; 
        // calculate the secondary sorption component concentrations
        ln_c_sec_sorp = _vec_lnK.segment( _J_mob, _J_sorp ) - Stoi_sorp * ln_c_basis; 
        lnK_min = _vec_lnK.tail(_J_min); 
    
        // fill in the secondary concentrations
        for (i=0; i < (size_t)ln_c_sec_mob.rows(); i++)
            c_second( i ) = std::exp( ln_c_sec_mob(i) );
        for (i=0; i < (size_t)ln_c_sec_sorp.rows(); i++)
            c_second( _I_sec_mob + i ) = std::exp( ln_c_sec_sorp(i) ); 
        c_second.tail( _I_sec_min )  = c_sec_min;
        // part 1), n_basis mass balance equations
        this->calc_tot_mass( c_basis, c_second, vec_cur_mass_balance ); 
        vec_residual.head( _I_basis ) = vec_tot_mass_constrain - vec_cur_mass_balance; 

        // part 2), n_react_min mineral reactions, 
        // AKA, the "complementary problem".
        for ( i=0; i < _J_min; i++ )
        {
            // attention, this is spectial for mineral reactions
            phi  = -1.0 * lnK_min(i) + Stoi_min.row(i) * ln_c_basis;
            res_tmp  = std::min( phi, c_sec_min(i) ); 
            vec_residual(_I_basis+i) = res_tmp; 
        }  // end of for
        */

    };  // end of function calc_residual
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
      * stoichiometric matrix from the input reactions
      */	
	LocalMatrix _matStoi_input;

    /**
      * vector of natural log K values
      */	
    LocalVector _vec_lnK; 
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
	size_t _I_mob, _I_sec_mob, _I_sec_sorp, _I_sec_min;

    /**
      * number of basis components
      */
	size_t _I_basis;

    /**
      * number of secondary components
      */
	size_t _I_second;

    /**
      * vector of system residual
      * its size is equal to the number of components
      */
    ogsChem::LocalVector _vec_res;

    /**
      * vector of staturation index
      * size of this vector is equal to the amount of minerals
      * if AI = 1, meaning the mineral is saturated and present
      * if AI = 0, meaning the mineral is under-saturated and not present
      */
    ogsChem::LocalVector _AI;

    /**
      * matrix of system Jacobi
      * its size is equal to the n_components by n_components
      */
    ogsChem::LocalMatrix _mat_Jacobi; 

    /**
      * count the components 
      */
    void countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
    {
	    _I_mob = 0; 
        _I_sec_mob = 0; 
	    _I_sec_sorp= 0; 
	    _I_sec_min = 0; 

	    BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
	    for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
	    {
		    switch ( it->second->getCompType() )
		    {
		    case ogsChem::AQ_PHASE_COMP: 
			    _I_mob++; 
			    break;
		    case ogsChem::SORPTION_COMP: 
			    _I_sec_sorp++;
			    break;
		    case ogsChem::MIN_PHASE_COMP: 
			    _I_sec_min++;
			    break;
		    default:
			    _I_sec_min++;
			    break; 
		    }
	    }

        _I = _I_mob + _I_sec_sorp + _I_sec_min;
    }; 

    /**
      * count how many mobile, sorption and mineral reactions
      */
    void countReactions(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
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
    };

    /**
      * construct stoichiometric matrix out of list of components and reactions
      */
    void buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
		           std::vector<ogsChem::chemReactionEq*>                & list_eq_reactions)
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
        _I_sec_mob = _I_mob - _I_basis; 
        // organize which components are basis
        // and which are secondary
        _matStoi = _matStoi_input.topRows( _I_basis ).transpose(); 

    #ifdef _DEBUG
	    // debugging--------------------------
	    std::cout << "Stoichiometric Matrix S: " << std::endl; 
	    std::cout << _matStoi << std::endl;
	    // end of debugging-------------------
    #endif
    };

    /**
      * reading in the logK values oaf each equilibrium reactions
      */
    void read_logK(std::vector<ogsChem::chemReactionEq*> & list_eq_reactions)
    {
        size_t i; 
        _vec_lnK = ogsChem::LocalVector::Zero(_J); 
        for ( i=0 ; i < list_eq_reactions.size(); i++ )
            _vec_lnK(i) = list_eq_reactions[i]->get_ln_K(); 
    };

};

}  // end of namespace

#endif