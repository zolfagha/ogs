
#pragma once

#include "NumLib/Coupling/CouplingStrucutreBuilder.h"
#include "NumLib/Coupling/Algorithm/TransientPartitionedAlgorithmFactory.h"

#include "TransientCoupledSystem.h"
#include "TransientMonolithicSystem.h"
#include "AsyncPartitionedSystem.h"


namespace NumLib
{

typedef class NumLib::TemplateCouplingStrucutreBuilder
	<
	ITransientCoupledSystem,
	TemplateTransientMonolithicSystem,
	AsyncPartitionedSystem,
	NumLib::TransientPartitionedAlgorithmFactory
	> TransientCoulplingStrucutreBuilder;

} //end
