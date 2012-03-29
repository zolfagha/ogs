
#pragma once

#include <vector>

namespace MathLib
{
class ICoupledSystem;
class ParameterSet;
struct ParameterProblemMappingTable;

/**
 * \brief
 */
class ITransientPartitionedAlgorithm
{
public:
	virtual ~ITransientPartitionedAlgorithm() {};

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<ICoupledSystem*> &subproblems, ParameterSet &vars_t_n, ParameterSet &vars_t_n1, ParameterProblemMappingTable &mapping) = 0;
};


}
