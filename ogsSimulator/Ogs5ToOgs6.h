
#pragma once

#include "FemIO/ogs5/Ogs5FemIO.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/PorousMedia.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "GeoLib/GEOObjects.h"

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertSolidProperty(const ogs5::CSolidProperties &msp, MaterialLib::Solid &solid);

void convertPorousMediumProperty(const ogs5::CMediumProperties &mmp, MaterialLib::PorousMedia &pm);

void convert(const ogs5::Ogs5FemData &ogs5fem, const GeoLib::GEOObjects &geo, const std::string &unique_geo_name, const DiscreteLib::DiscreteSystem &dis);

};
} //end
