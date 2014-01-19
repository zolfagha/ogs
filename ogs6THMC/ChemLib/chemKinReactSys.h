/**
* Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file chemKinReactSys.h
*
* Created on 2014-01-19 by Haibing Shao
*/
#ifndef CHEM_KIN_REACT_SYS_H
#define CHEM_KIN_REACT_SYS_H

#include "chemconst.h"
#include "chemcomp.h"
#include "chemReactionKin.h"
#include "BaseLib/OrderedMap.h"
#include "chemActivityModelUnity.h"

namespace ogsChem
{

class chemKinReactSys
{
public:
    /**
    * constructor of the class
    */
    chemKinReactSys(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp,
                    std::vector<ogsChem::chemReactionKin*>               & list_kin_reactions,
                    ogsChem::chemActivityModelAbstract                  *a)
                    : _list_kin_reactions(list_kin_reactions), 
                      _I(map_chemComp.size()), _J_kin(list_kin_reactions.size()),
                      _activity_model(a)
    {
        if (_I > 0 && _J_kin > 0)
            isInitialized = true; 
    
    };

    /**
      * destructor of the class
      */
    ~chemKinReactSys(void)
    {};


    /**
      * solve the kinetic reaction system
      * using the StepperBulischStoer ODE solver. 
      * input values are: 
      * the vector of concentration values
      * the convergence tolerance,
      * the maximum number of iterations.
      * and the status of convergence
      */
    void solve_KinSys(LocalVector & vec_conc,
                      size_t & result,
                      size_t & node_idx,
                      double iter_tol,
                      double rel_tol,
                      double dt, 
                      double max_iter)
    {
        
    }; 

private:

    /**
      * private flag indicating initialization
      */
    bool isInitialized;

    /**
      * a list of all kinetic reactions
      */
    std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions;

    /**
      * pointer to the activity model
      */
    ogsChem::chemActivityModelAbstract *_activity_model;

    /**
      * _I is the number of components and _J_kin is the number of kinetic reactions
      */
    size_t _I, _J_kin;

};

}

#endif