
#pragma once

#include <vector>

#include "CouplingSolution.h"

namespace NumLib
{
typedef std::pair<ICouplingProblem*,size_t> PairSysVarId;

class IPartitionedAlgorithm
{
public:
    virtual int solve(std::vector<ICouplingProblem*> &list_children, std::vector<std::pair<std::string, PairSysVarId>> &map_sharedVar2sysVar, SharedVariables* vars) = 0;
};

class IterativePartitionedMethod : public IPartitionedAlgorithm
{
public:
    IterativePartitionedMethod() : _max_itr(200)
    {
    }

    size_t getIterationCounts() const {return _itr_count;};
    void setMaximumIterationCounts(size_t n) {_max_itr = n;};
    size_t getMaximumIterationCounts() const {return _max_itr;};

protected:
    void setIterationCounts(size_t n) {_itr_count = n;};

private:
    size_t _max_itr;
    size_t _itr_count;
};


/**
 * \brief Block Jacobi method
 */
class BlockJacobiMethod : public IterativePartitionedMethod
{
public:
    int solve(std::vector<ICouplingProblem*> &list_children, std::vector<std::pair<std::string, PairSysVarId>> &map_sharedVar2sysVar, SharedVariables* vars)
    {
        size_t max_itr = getMaximumIterationCounts();

        for (size_t i_itr=0; i_itr<max_itr; i_itr++) {
            // set previous state
            for (size_t i=0; i<map_sharedVar2sysVar.size(); i++) {
                const std::string &name = map_sharedVar2sysVar[i].first;
                ICouplingProblem *solution = map_sharedVar2sysVar[i].second.first;
                size_t internalId = map_sharedVar2sysVar[i].second.second;
                vars->setVariable(name, const_cast<Variable*>(solution->get(internalId)));
            }
            // compute each
            for (size_t i=0; i<list_children.size(); i++) {
                ICouplingProblem *solution = list_children[i];
                solution->solve();
            }
            // check convergence
            if (isConverged(vars)) {
                setIterationCounts(i_itr);
                break;
            }
        }

        return 0;
    }

    bool isConverged(SharedVariables* vars)
    {
        return true;
    }
};


class BlockGaussSeidelMethod : public IterativePartitionedMethod
{

};


class PartitionedStaggeredMethod : public IPartitionedAlgorithm
{

};

class ParallelStaggeredMethod : public PartitionedStaggeredMethod
{

};

class SerialStaggeredMethod : public PartitionedStaggeredMethod
{

};
}
