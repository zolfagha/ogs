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
#include "UserModules/FemReactOPS/Local_ODE_KinReact.h"
#include "UserModules/FemReactGIAReduct/StepperBulischStoer.h"

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
        double t0 = 0.0;

        _p_local_ODE = new Local_ODE_KinReact(list_kin_reactions);

        _vec_Comp_Conc = ogsChem::LocalVector::Zero(_I); 

        _vec_xi_kin_rate = ogsChem::LocalVector::Zero(_J_kin); 

        _sbs = new MathLib::StepperBulischStoer<Local_ODE_KinReact>(_vec_Comp_Conc, _vec_xi_kin_rate,
                                                                    t0, 1.0e-12, 1.0e-12, true);

        if (_I > 0 && _J_kin > 0 && _p_local_ODE && _sbs)
            isInitialized = true; 

    };

    /**
      * destructor of the class
      */
    ~chemKinReactSys(void)
    {
        BaseLib::releaseObject(_p_local_ODE);
        BaseLib::releaseObject(_sbs);
    };


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
        // getting the initial rate evaluation. 
        if (_vec_xi_kin_rate.rows() > 0)
            _vec_xi_kin_rate = (*_p_local_ODE)(dt, vec_conc);

        //_sbs->set_y(Xi_Kin_bar);
        _sbs->set_y(vec_conc);
        _sbs->set_dydx(_vec_xi_kin_rate);

        // solve the local ODE problem for xi kin bar
        _sbs->step(dt, _p_local_ODE);
        //vec_Xi_Kin_bar_new = _sbs->get_y();
        vec_conc = _sbs->get_y();
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

    /**
      * pointer to the Local_ODE_KinReact
      */
    Local_ODE_KinReact* _p_local_ODE; 

    /**
      * pointer to the ODE solver StepperBulischStoer
      */
    MathLib::StepperBulischStoer<Local_ODE_KinReact>* _sbs;

    ogsChem::LocalVector _vec_Comp_Conc; 

    ogsChem::LocalVector _vec_xi_kin_rate; 

};

}

#endif