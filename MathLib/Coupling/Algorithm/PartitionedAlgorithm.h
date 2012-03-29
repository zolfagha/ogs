
#pragma once

#include <vector>

namespace MathLib
{
class ICoupledSystem;
class ParameterSet;
struct ParameterProblemMappingTable;

/**
 * \brief Interface class for partitioned algorithm
 */
class IPartitionedAlgorithm
{
public:

	virtual ~IPartitionedAlgorithm() {};

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(const std::vector<ICoupledSystem*> &subproblems, ParameterSet &vars, const ParameterProblemMappingTable &mapping) = 0;
};

}
