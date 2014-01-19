/**
* Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.com/LICENSE.txt
*
*
* \file Local_ODE_Xi_immob_GIA.h
*
* Created on 2014-Jan-19 by Haibing Shao
*/

#ifndef LOCAL_ODE_KIN_REACT_H
#define LOCAL_ODE_KIN_REACT_H

#include "ChemLib/chemcomp.h"
#include "ChemLib/chemReactionKin.h""

class Local_ODE_KinReact
{
public:

    /**
    * constructor
    */
    Local_ODE_KinReact(std::vector<ogsChem::chemReactionKin*> & list_kin_reactions)
                       : _list_kin_reactions(list_kin_reactions)
    {
        _J_kin = _list_kin_reactions.size();

        _vec_dxi_immob_dt_new = ogsChem::LocalVector::Zero(_J_kin); 
    }

    /**
      * destructor
      */
    ~Local_ODE_KinReact()
    {
    }

    /**
      * evaluate the change of xi_immob over time
      */
    ogsChem::LocalVector operator() (double dt, MathLib::LocalVector & vec_Comp_Conc)
    {
        for (std::size_t i = 0; i < _J_kin; i++)
        {
            _list_kin_reactions[i]->calcReactionRate(vec_Comp_Conc);
            _vec_dxi_immob_dt_new(i) = _list_kin_reactions[i]->getRate(); 
        }
        return _vec_dxi_immob_dt_new;
    }

private:

    /**
      * vector storing the kinetic reaction rates
      */
    ogsChem::LocalVector _vec_dxi_immob_dt_new; 

    /**
      * address of kinetic reactions
      */
    std::vector<ogsChem::chemReactionKin*> & _list_kin_reactions; 

    /**
      * number of kinetic reactions
      */
    std::size_t _J_kin; 
};

#endif
