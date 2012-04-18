
#pragma once

#include "MathLib/Coupling/CouplingStrucutreBuilder.h"
#include "MathLib/Coupling/Algorithm/TransientPartitionedAlgorithmFactory.h"

#include "TransientCoupledSystem.h"
#include "TransientMonolithicSystem.h"
#include "AsyncPartitionedSystem.h"


namespace NumLib
{

typedef class MathLib::TemplateCouplingStrucutreBuilder
	<
	ITransientCoupledSystem,
	TemplateTransientMonolithicSystem,
	AsyncPartitionedSystem,
	MathLib::TransientPartitionedAlgorithmFactory
	> TransientCoulplingStrucutreBuilder;

} //end
