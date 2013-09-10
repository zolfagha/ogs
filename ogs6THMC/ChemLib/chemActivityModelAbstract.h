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
     * set current time
     * @param molarity
	 * @param activity
     */
    virtual void calc_activity(MathLib::LocalVector molarity, MathLib::LocalVector activity) = 0;

protected:
    ogsChem::ACTIVITY_MODEL  _activity_model; 
};

}

#endif