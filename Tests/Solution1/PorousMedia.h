
#pragma once

#include "Base/CodingTools.h"
#include "MathLib/Function/IFunction.h"

struct PorousMedia
{
	MathLib::TemplateFunction<double*, double>* K;
	MathLib::TemplateFunction<double*, double>* porosity;

	PorousMedia()
	{
		Base::zeroObject(K, porosity);
	}
	~PorousMedia()
	{
		Base::releaseObject(K, porosity);
	}
};
