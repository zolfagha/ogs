/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Compound.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"

namespace MaterialLib
{

struct Compound
{
    std::string name;
    bool is_mobile;
    NumLib::ITXFunction* molecular_diffusion;

    Compound()
    {
        BaseLib::zeroObject(molecular_diffusion);
    }
    ~Compound()
    {
        BaseLib::releaseObject(molecular_diffusion);
    }

};

} //end

