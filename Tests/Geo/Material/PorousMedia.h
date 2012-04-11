
#pragma once

#include "Base/CodingTools.h"
#include "MathLib/Function/IFunction.h"

namespace Geo
{

struct PorousMedia
{
	MathLib::SpatialFunctionScalar* hydraulic_conductivity;
	MathLib::SpatialFunctionScalar* porosity;
	MathLib::SpatialFunctionScalar* storage;

	PorousMedia()
	{
		Base::zeroObject(hydraulic_conductivity, porosity, storage);
	}
	~PorousMedia()
	{
		Base::releaseObject(hydraulic_conductivity, porosity, storage);
	}
};

} //end
 