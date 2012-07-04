
#pragma once

#include "BaseLib/CodingTools.h"

#include "Solid.h"

namespace NumLib
{
class ITXFunction;
}

namespace MaterialLib
{

struct PorousMedia
{
	NumLib::ITXFunction* hydraulic_conductivity;
	NumLib::ITXFunction* porosity;
	NumLib::ITXFunction* storage;
	Solid* solidphase;

	PorousMedia()
	{
		BaseLib::zeroObject(hydraulic_conductivity, porosity, storage, solidphase);
	}
	~PorousMedia()
	{
		BaseLib::releaseObject(hydraulic_conductivity, porosity, storage);
	}
};

} //end
 
