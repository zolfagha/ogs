
#pragma once

#include "MathLib/IOSystem/ParameterSet.h"

namespace MathLib
{

class IConvergenceCheck
{
public:
	virtual ~IConvergenceCheck() {};

	virtual bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff) = 0;
};

} //end
