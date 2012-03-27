
#pragma once

#include "IterativePartitionedMethod.h"

namespace NumLib
{

/**
 * \brief Block Gauss-Seidel iterative partitioned method
 */
template <class T_CONVERGENCE_CHECK>
class BlockGaussSeidelMethod : public AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>
{
public:
	typedef MathLib::IFunction Variable;
	typedef ICoupledSystem MyCoupledSystem;
	typedef VariableContainer MyNamedVariableContainer;
    typedef VariableMappingTable MyVariableMappingTable;

    ///
    BlockGaussSeidelMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>(epsilon, max_count)
    {
    }
    ///
    virtual ~BlockGaussSeidelMethod() {};

    ///
    void doPostAfterSolve( MyCoupledSystem & solution, MyNamedVariableContainer& vars, MyVariableMappingTable &mapping );
};

template <class T_CONVERGENCE_CHECK>
void BlockGaussSeidelMethod<T_CONVERGENCE_CHECK>::doPostAfterSolve( MyCoupledSystem & solution, MyNamedVariableContainer& vars, MyVariableMappingTable &mapping )
{
    // update shared variables
    const size_t n_vars = vars.size();
    for (size_t i=0; i<n_vars; i++) {
        typename MyVariableMappingTable::PairSysVarId &v = mapping._map_paraId2subproblem[i];
        const MyCoupledSystem *work_solution = v.first;
        if (work_solution==&solution) {
            vars.set(i, *const_cast<Variable*>(solution.getParameter(v.second)));
        }
    }
}


}
