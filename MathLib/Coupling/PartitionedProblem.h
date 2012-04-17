
#pragma once

#include <vector>
#include <iostream>

#include "MathLib/Parameter/ParameterSet.h"
#include "MathLib/Parameter/NamedIOSystem.h"

#include "ICoupledProblem.h"
#include "AbstractPartitionedProblem.h"
#include "ParameterProblemMappingTable.h"
#include "Algorithm/PartitionedAlgorithm.h"

namespace MathLib
{

/**
 * \brief Partitioned problem
 */
class PartitionedProblem : public AbstractPartitionedProblem<ICoupledSystem>
{
public:
	typedef size_t InternalID;
	typedef size_t ExternalKey;

    /// 
    PartitionedProblem() : _algorithm(0)
    {
    }

    ///
    virtual ~PartitionedProblem() {};

    ///
    virtual void setAlgorithm(IPartitionedAlgorithm &algo)
    {
    	_algorithm = &algo;
    }

    /// solve this system
    virtual int solve();

private:
    IPartitionedAlgorithm *_algorithm;
};



}
