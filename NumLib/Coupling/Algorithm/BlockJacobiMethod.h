
#pragma once

#include "AbstractIterativePartitionedMethod.h"

namespace NumLib
{

/**
 * \brief Block Jacobi iterative partitioned method
 */
class BlockJacobiMethod : public AbstractIterativePartitionedMethod
{
public:
	BlockJacobiMethod() {};

    BlockJacobiMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(epsilon, max_count)
    {
    }

    BlockJacobiMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(checker, epsilon, max_count)
    {
    }

    virtual ~BlockJacobiMethod() {};

protected:
    void doPostAfterSolveAll( UnnamedParameterSet &parameter_table, const ParameterProblemMappingTable &mapping ) //;
//    void BlockJacobiMethod::doPostAfterSolveAll( UnnamedParameterSet &parameter_table, const ParameterProblemMappingTable &mapping )
    {
        // set current state to shared variables
    	AbstractIterativePartitionedMethod::updateParameterTable(mapping, isFixed(), parameter_table);
    }

    bool isFixed() const {return false;};
};


} //end