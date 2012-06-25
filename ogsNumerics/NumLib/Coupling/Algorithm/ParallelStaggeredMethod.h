
#pragma once

#include "AbstractPartitionedStaggeredMethod.h"
#include "BlockJacobiMethod.h"

namespace NumLib
{

class ParallelStaggeredMethod : public AbstractPartitionedStaggeredMethod<NumLib::BlockJacobiMethod >
{
public:
    ParallelStaggeredMethod() {};
    ParallelStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockJacobiMethod >(epsilon, max_count) {};
    ParallelStaggeredMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockJacobiMethod >(checker, epsilon, max_count) {};
    virtual ~ParallelStaggeredMethod() {};
};

}
