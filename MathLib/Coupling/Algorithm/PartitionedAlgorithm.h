
#pragma once

#include <vector>
#include "MathLib/Coupling/ICoupledProblem.h"

namespace MathLib
{
class UnnamedParameterSet;
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
    virtual int solve(const std::vector<ICoupledSystem*> &subproblems, UnnamedParameterSet &vars, const ParameterProblemMappingTable &mapping) = 0;
};

}
