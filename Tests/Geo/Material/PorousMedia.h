
#pragma once

#include "Base/CodingTools.h"
#include "NumLib/Function/IFunction.h"

#include "Solid.h"

namespace Geo
{

struct PorousMedia
{
	NumLib::SpatialFunctionScalar* hydraulic_conductivity;
	NumLib::SpatialFunctionScalar* porosity;
	NumLib::SpatialFunctionScalar* storage;
	Solid* solidphase;

	PorousMedia()
	{
		Base::zeroObject(hydraulic_conductivity, porosity, storage, solidphase);
	}
	~PorousMedia()
	{
		Base::releaseObject(hydraulic_conductivity, porosity, storage, solidphase);
	}
};

} //end
 
