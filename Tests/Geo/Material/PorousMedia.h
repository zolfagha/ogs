
#pragma once

#include "Base/CodingTools.h"
#include "MathLib/Function/IFunction.h"

#include "Solid.h"

namespace Geo
{

struct PorousMedia
{
	MathLib::SpatialFunctionScalar* hydraulic_conductivity;
	MathLib::SpatialFunctionScalar* porosity;
	MathLib::SpatialFunctionScalar* storage;
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
 
