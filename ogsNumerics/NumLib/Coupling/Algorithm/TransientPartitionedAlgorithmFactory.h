
#pragma once

#include <string>
#include "TransientPartitionedAlgorithm.h"
#include "IConvergenceCheck.h"

namespace NumLib
{

class TransientPartitionedAlgorithmFactory
{
public:
    static ITransientPartitionedAlgorithm* create(const std::string &name, IConvergenceCheck *checker, size_t max_itr, double epsilon)
    {
        if (name.compare("Parallel")==0) {
            return new ParallelStaggeredMethod(*checker, epsilon, max_itr);
        } else if (name.compare("Serial")==0) {
            return new SerialStaggeredMethod(*checker, epsilon, max_itr);
        }
        return 0;
    };
};


} //end
