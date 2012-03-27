
#pragma once

#include "IterativePartitionedMethod.h"

namespace NumLib
{

/**
 * \brief Block Jacobi iterative partitioned method
 */
template <class T_CONVERGENCE_CHECK>
class BlockJacobiMethod : public AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>
{
public:
	typedef MathLib::IFunction Variable;
	typedef ICoupledSystem MyCoupledSystem;
	typedef VariableContainer MyNamedVariableContainer;
    typedef VariableMappingTable MyVariableMappingTable;
    BlockJacobiMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>(epsilon, max_count)
    {
    }
    virtual ~BlockJacobiMethod() {};

    void doPostAfterSolveAll( MyNamedVariableContainer &vars, MyVariableMappingTable &mapping );
};

template <class T_CONVERGENCE_CHECK>
void BlockJacobiMethod<T_CONVERGENCE_CHECK>::doPostAfterSolveAll( MyNamedVariableContainer &vars, MyVariableMappingTable &mapping )
{
    // set current state to shared variables
    const size_t n_vars = vars.size();
    for (size_t i=0; i<n_vars; i++) {
    	typename MyVariableMappingTable::PairSysVarId &v = mapping._map_paraId2subproblem[i];
        const MyCoupledSystem *solution = v.first;
        if (solution!=0)
            vars.set(i, *const_cast<Variable*>(solution->getParameter(v.second)));
    }
}

}
