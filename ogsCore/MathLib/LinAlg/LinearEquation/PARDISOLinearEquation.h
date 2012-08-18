/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PARDISOInterface.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include "AbstractCRSLinearEquation.h"

namespace MathLib
{

class PARDISO_Solver : public AbstractCRSLinearEquation<signed>
{

};

}
