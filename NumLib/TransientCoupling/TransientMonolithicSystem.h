
#pragma once

#include <cassert>
#include <vector>

#include "Base/CodingTools.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

//template <size_t N_IN, size_t N_OUT>
//class TemplateTransientMonolithicSystem : public NumLib::TemplateMonolithicSystem<ITransientCoupledSystem, N_IN, N_OUT>
class TemplateTransientMonolithicSystem : public NumLib::AbstractMonolithicSystem<ITransientCoupledSystem>
{

};


}
