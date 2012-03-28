
#pragma once

#include "IterativePartitionedMethod.h"

namespace MathLib
{

/**
 * \brief Block Gauss-Seidel iterative partitioned method
 */
template <class T_CONVERGENCE_CHECK>
class BlockGaussSeidelMethod : public AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>
{
public:

    ///
    BlockGaussSeidelMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>(epsilon, max_count)
    {
    }
    ///
    virtual ~BlockGaussSeidelMethod() {};

    ///
    void doPostAfterSolve( ICoupledSystem & problem, VariableContainer& parameter_table, const VariableMappingTable &mapping );
};

template <class T_CONVERGENCE_CHECK>
void BlockGaussSeidelMethod<T_CONVERGENCE_CHECK>::doPostAfterSolve( ICoupledSystem & problem, VariableContainer& parameter_table, const VariableMappingTable &mapping )
{
    // update shared variables
	AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>::updateParameterTable(problem, mapping, parameter_table);
}


}
