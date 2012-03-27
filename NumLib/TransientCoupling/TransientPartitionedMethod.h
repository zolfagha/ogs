
#pragma once

#include "NumLib/Coupling/VariableMappingTable.h"
#include "NumLib/Coupling/Algorithm/BlockGaussSeidelMethod.h"
#include "NumLib/Coupling/Algorithm/BlockJacobiMethod.h"
#include "TransientPartitionedAlgorithm.h"

namespace NumLib
{
/**
 * \brief 
 */
template <class T_PARTITIONED>
class AbstractPartitionedStaggeredMethod : public ITransientPartitionedAlgorithm
{
public:
	typedef VariableContainer MyNamedVariableContainer;
	typedef MathLib::IFunction Variable;
	typedef ICoupledSystem MyCoupledSystem;
    typedef VariableMappingTable MyVariableMappingTable;

    AbstractPartitionedStaggeredMethod() : _part_method() {};
    AbstractPartitionedStaggeredMethod(double epsilon, size_t max_count) : _part_method(epsilon, max_count) {};
    virtual ~AbstractPartitionedStaggeredMethod() {};

    /// solve
    int solve(std::vector<MyCoupledSystem*> &subproblems, MyNamedVariableContainer &vars_t_n, MyNamedVariableContainer &vars_t_n1, MyVariableMappingTable &mapping)
    {
        // initialize variables at t_n1 from t_n
        const size_t n_para = vars_t_n.size();
        for (size_t i=0; i<n_para; i++) {
            typename MyVariableMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            MyCoupledSystem *solution = sysvarid.first;
            //size_t internalId = sysvarid.second;
            if (solution!=0)
                vars_t_n1.set(i, *vars_t_n.get(i));
        }

        _part_method.solve(subproblems, vars_t_n1, mapping);

        return 0;
    }

private:
    T_PARTITIONED _part_method;
};

template <class T_CONVERGENCE_CHECK>
class ParallelStaggeredMethod : public AbstractPartitionedStaggeredMethod<BlockJacobiMethod<T_CONVERGENCE_CHECK> >
{
public:
	ParallelStaggeredMethod() {};
	ParallelStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockJacobiMethod<T_CONVERGENCE_CHECK> >(epsilon, max_count) {};
    virtual ~ParallelStaggeredMethod() {};
};

template <class T_CONVERGENCE_CHECK>
class SerialStaggeredMethod : public AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod<T_CONVERGENCE_CHECK> >
{
public:
	SerialStaggeredMethod() {};
	SerialStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod<T_CONVERGENCE_CHECK> >(epsilon, max_count) {};
    virtual ~SerialStaggeredMethod() {};
};

//typedef AbstractPartitionedStaggeredMethod<BlockJacobiMethod> ParallelStaggeredMethod;
//typedef AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod> SerialStaggeredMethod;

}
