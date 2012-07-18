
#pragma once

#include "BaseLib/CodingTools.h"


namespace NumLib
{
class ITXFunction;
}

namespace MaterialLib
{

struct PorousMedia
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* permeability;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;

    PorousMedia()
    {
        BaseLib::zeroObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage
                );
    }

    ~PorousMedia()
    {
        BaseLib::releaseObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage
                );
    }
};

} //end
 
