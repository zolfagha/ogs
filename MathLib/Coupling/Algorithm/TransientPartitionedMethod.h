
#pragma once

#include "MathLib/Coupling/ParameterTable.h"
#include "BlockGaussSeidelMethod.h"
#include "BlockJacobiMethod.h"
#include "TransientPartitionedAlgorithm.h"

namespace MathLib
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
    virtual ~AbstractPartitionedStaggeredMethod() {};

    /// solve
    int solve(std::vector<ICoupledSystem*> &subproblems, ParameterSet &parameters_t_n, ParameterSet &parameters_t_n1, ParameterProblemMappingTable &mapping)
    {
        initializeParameters(parameters_t_n, subproblems, mapping, parameters_t_n1);

        _part_method.solve(subproblems, parameters_t_n1, mapping);

        return 0;
    }

private:
    void initializeParameters(const ParameterSet &parameters_t_n, const std::vector<ICoupledSystem*> &subproblems, const ParameterProblemMappingTable &mapping, ParameterSet &parameters_t_n1)
    {
        const size_t n_para = parameters_t_n.size();
        for (size_t i=0; i<n_para; i++) {
            const ParameterProblemMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            const ICoupledSystem *solution = sysvarid.first;
            //size_t internalId = sysvarid.second;
            if (solution!=0)
                parameters_t_n1.set(i, *parameters_t_n.get(i));
        }
    }
private:
    T_PARTITIONED _part_method;
};

template <class T_CONVERGENCE_CHECK>
class ParallelStaggeredMethod : public AbstractPartitionedStaggeredMethod<MathLib::BlockJacobiMethod<T_CONVERGENCE_CHECK> >
{
public:
    ParallelStaggeredMethod() {};
    ParallelStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockJacobiMethod<T_CONVERGENCE_CHECK> >(epsilon, max_count) {};
    virtual ~ParallelStaggeredMethod() {};
};

template <class T_CONVERGENCE_CHECK>
class SerialStaggeredMethod : public AbstractPartitionedStaggeredMethod<MathLib::BlockGaussSeidelMethod<T_CONVERGENCE_CHECK> >
{
public:
    SerialStaggeredMethod() {};
    SerialStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod<T_CONVERGENCE_CHECK> >(epsilon, max_count) {};
    virtual ~SerialStaggeredMethod() {};
};

}
