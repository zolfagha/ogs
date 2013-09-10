/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemActivityModelAbstract.h
 *
 * Created on 2013-09-10 by Haibing Shao
 */

#ifndef CHEM_ACTIVITY_MODEL_ABSTRACT_H
#define CHEM_ACTIVITY_MODEL_ABSTRACT_H

#include "chemconst.h"

namespace ogsChem
{

/**
 * \brief Common interface of Activity Model
 */
class chemActivityModelAbstract
{
public:
    // destructor
    virtual ~chemActivityModelAbstract() {};

    /**
     * calculate the activity based on molarity
     * @param log_molarity
	 * @param log_activity
     */
    virtual void calc_activity(MathLib::LocalVector & log_molarity, 
                               MathLib::LocalVector & log_activity, 
                               const double Tc = 25.0) = 0;

protected:
    const ogsChem::ACTIVITY_MODEL  _activity_model; 
};

} // end of namespace

#endif