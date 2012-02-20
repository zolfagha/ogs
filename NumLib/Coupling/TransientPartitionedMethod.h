
#pragma once

#include "IterativePartitionedMethod.h"

#include "VariableMappingTable.h"

namespace NumLib
{
/**
 * \brief 
 */
template <class T_PARTITIONED>
class AbstractPartitionedStaggeredMethod : public ITransientPartitionedAlgorithm
{
public:
    AbstractPartitionedStaggeredMethod() : _part_method() {};
    AbstractPartitionedStaggeredMethod(double epsilon, size_t max_count) : _part_method(epsilon, max_count) {};

    /// solve
    int solve(std::vector<ICoupledSystem*> &subproblems, NamedVariableContainer &vars_t_n, NamedVariableContainer &vars_t_n1, VariableMappingTable &mapping)
    {
        // initialize variables at t_n1 from t_n
        const size_t n_para = vars_t_n.size();
        for (size_t i=0; i<n_para; i++) {
            VariableMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            ICoupledSystem *solution = sysvarid.first;
            size_t internalId = sysvarid.second;
            if (solution!=0)
                vars_t_n1.set(i, *vars_t_n.get(i));
        }

        _part_method.solve(subproblems, vars_t_n1, mapping);

        return 0;
    }

private:
    T_PARTITIONED _part_method;
};

typedef AbstractPartitionedStaggeredMethod<BlockJacobiMethod> ParallelStaggeredMethod;
typedef AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod> SerialStaggeredMethod;

}
