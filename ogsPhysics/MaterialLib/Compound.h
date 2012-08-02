/*
 * Compound.h
 *
 *  Created on: 10.04.2012
 *      Author: watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

namespace NumLib
{
class ITXFunction;
}

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

