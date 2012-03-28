
#pragma once

#include "IterativePartitionedMethod.h"

namespace MathLib
{

/**
 * \brief Block Jacobi iterative partitioned method
 */
template <class T_CONVERGENCE_CHECK>
class BlockJacobiMethod : public AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>
{
public:
    BlockJacobiMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>(epsilon, max_count)
    {
    }
    virtual ~BlockJacobiMethod() {};

    void doPostAfterSolveAll( ParameterTable &parameter_table, const ParameterProblemMappingTable &mapping );
};

template <class T_CONVERGENCE_CHECK>
void BlockJacobiMethod<T_CONVERGENCE_CHECK>::doPostAfterSolveAll( ParameterTable &parameter_table, const ParameterProblemMappingTable &mapping )
{
    // set current state to shared variables
	AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>::updateParameterTable(mapping, parameter_table);
}

} //end
