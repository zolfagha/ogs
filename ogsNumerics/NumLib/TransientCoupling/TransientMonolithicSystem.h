
#pragma once

#include <cassert>
#include <vector>

#include "BaseLib/CodingTools.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

class AbstractTransientMonolithicSystem : public NumLib::AbstractMonolithicSystem<ITransientCoupledSystem>
{
public:
	virtual ~AbstractTransientMonolithicSystem() {};
};


}
