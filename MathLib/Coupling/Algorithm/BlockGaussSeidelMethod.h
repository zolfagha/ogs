
#pragma once

#include "AbstractIterativePartitionedMethod.h"

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

protected:
    bool isFixed() const {return true;};
    ///
    void doPostAfterSolve( ICoupledSystem & problem, ParameterSet& parameter_table, const ParameterProblemMappingTable &mapping );
};

template <class T_CONVERGENCE_CHECK>
void BlockGaussSeidelMethod<T_CONVERGENCE_CHECK>::doPostAfterSolve( ICoupledSystem & problem, ParameterSet& parameter_table, const ParameterProblemMappingTable &mapping )
{
    // update shared variables
	AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>::updateParameterTable(problem, mapping, isFixed(), parameter_table);
}


}
