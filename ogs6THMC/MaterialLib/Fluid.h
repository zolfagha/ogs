/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Fluid.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"

namespace MaterialLib
{

struct Fluid
{
    NumLib::ITXFunction* density;
    NumLib::ITXFunction* dynamic_viscosity;
    NumLib::ITXFunction* specific_heat;
    NumLib::ITXFunction* thermal_conductivity;
	NumLib::ITXFunction* drho_dp;

    Fluid()
    {
        BaseLib::zeroObject(
                density,
                dynamic_viscosity,
                specific_heat,
                thermal_conductivity,
				drho_dp
                );
    }
    ~Fluid()
    {
        BaseLib::releaseObject(
                density,
                dynamic_viscosity,
                specific_heat,
                thermal_conductivity,
				drho_dp
                );
    }
};

} //end
