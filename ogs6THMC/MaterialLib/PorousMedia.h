/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "IMedium.h"

namespace NumLib
{
class ITXFunction;
}

namespace MaterialLib
{

struct PorousMedia : public IMedium
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* permeability;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;
    NumLib::ITXFunction* geo_area;
    NumLib::ITXFunction* dispersivity_long;
    NumLib::ITXFunction* dispersivity_trans; 

    PorousMedia()
    {
        BaseLib::zeroObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area
                );
    }

    virtual ~PorousMedia()
    {
        BaseLib::releaseObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area
                );
    }

    virtual MediumType getMediumType() const {return MediumType::PorousMedium;};
};

} //end
 
