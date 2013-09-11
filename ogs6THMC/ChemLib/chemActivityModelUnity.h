/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemActivityModelAqUnity.h
 *
 * Created on 2013-09-10 by Haibing Shao
 */

#ifndef CHEM_ACTIVITY_MODEL_UNITY_H
#define CHEM_ACTIVITY_MODEL_UNITY_H

#include "chemActivityModelAbstract.h"

namespace ogsChem
{

/**
 * \brief Common interface of Activity Model
 */
class chemActivityModelUnity : chemActivityModelAbstract
{
public:
    // destructor
    chemActivityModelUnity() 
        : chemActivityModelAbstract(ogsChem::ACT_MOD_UNITY)
    {};

    /**
     * calculate the activity based on molarity
     * @param log_molarity
	 * @param log_activity
     */
    void calc_activity(MathLib::LocalVector & log_molarity, 
                       MathLib::LocalVector & log_activity, 
                       const double Tc = 25.0)
    {
        log_activity = log_molarity; 
    }

};

} // end of namespace

#endif