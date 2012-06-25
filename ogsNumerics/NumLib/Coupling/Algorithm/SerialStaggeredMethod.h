
#pragma once

#include "AbstractPartitionedStaggeredMethod.h"
#include "BlockGaussSeidelMethod.h"

namespace NumLib
{

class SerialStaggeredMethod : public AbstractPartitionedStaggeredMethod<NumLib::BlockGaussSeidelMethod >
{
public:
    SerialStaggeredMethod() {};
    SerialStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod >(epsilon, max_count) {};
    SerialStaggeredMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod >(checker, epsilon, max_count) {};
    virtual ~SerialStaggeredMethod() {};
};

}
