
#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"

namespace MaterialLib
{

struct Fluid
{
    NumLib::ITXFunction* density;
    NumLib::ITXFunction* dynamic_viscosity;

    Fluid()
    {
        BaseLib::zeroObject(
                density,
                dynamic_viscosity
                );
    }
    ~Fluid()
    {
        BaseLib::releaseObject(
                density,
                dynamic_viscosity
                );
    }
};

} //end
