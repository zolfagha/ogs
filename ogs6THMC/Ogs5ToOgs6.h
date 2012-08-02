
#pragma once

#include "FemIO/ogs5/Ogs5FemIO.h"
#include "MaterialLib/Solid.h"
#include "MaterialLib/PorousMedia.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "GeoLib/GEOObjects.h"
#include "Ogs6FemData.h"

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertSolidProperty(const ogs5::CSolidProperties &msp, MaterialLib::Solid &solid);

void convertPorousMediumProperty(const ogs5::CMediumProperties &mmp, MaterialLib::PorousMedia &pm);

void convert(const ogs5::Ogs5FemData &ogs5fem, Ogs6FemData &ogs6fem, BaseLib::Options &option);
};
} //end
