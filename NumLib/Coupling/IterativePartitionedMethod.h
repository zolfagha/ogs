
#pragma once

#include "PartitionedAlgorithm.h"

namespace NumLib
{

/**
 * \brief Abstract class for iterative partitioned methods
 */
class AbstractIterativePartitionedMethod : public IPartitionedAlgorithm
{
public:
    AbstractIterativePartitionedMethod() : _max_itr(100), _epsilon(1e-3)
    {
    }

    AbstractIterativePartitionedMethod(double epsilon, size_t max_count) : _max_itr(max_count), _epsilon(epsilon)
    {
    }

    /// get max. iteration
    size_t getMaximumIterationCounts() const {return _max_itr;};
    /// get epsilon
    double getEpsilon() const {return _epsilon;};
    /// get iteration count
    size_t getIterationCounts() const {return _itr_count;};
    /// solve
    int solve(std::vector<ICoupledProblem*> &subproblems, NamedVariableContainer &vars, VariableMappingTable &mapping);

protected:
    /// check if solution is converged
    bool isConverged(NamedVariableContainer& vars_prev, NamedVariableContainer& vars_current, double &v_diff);

    virtual void doPostAfterSolve( ICoupledProblem& solution, NamedVariableContainer& vars, VariableMappingTable &mapping )  {}

    virtual void doPostAfterSolveAll( NamedVariableContainer &vars, VariableMappingTable &mapping ) {}

private:
    size_t _max_itr;
    size_t _itr_count;
    double _epsilon;
};


/**
 * \brief Block Jacobi iterative partitioned method
 */
class BlockJacobiMethod : public AbstractIterativePartitionedMethod
{
public:
    BlockJacobiMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(epsilon, max_count)
    {
    }

    void doPostAfterSolveAll( NamedVariableContainer &vars, VariableMappingTable &mapping );
};


/**
 * \brief Block Gauss-Seidel iterative partitioned method
 */
class BlockGaussSeidelMethod : public AbstractIterativePartitionedMethod
{
public:
    BlockGaussSeidelMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(epsilon, max_count)
    {
    }

    void doPostAfterSolve( ICoupledProblem & solution, NamedVariableContainer& vars, VariableMappingTable &mapping );
};

}
