
#pragma once

#include <vector>
#include "NumLib/Coupling/ICoupledProblem.h"
#include "IConvergenceCheck.h"

namespace NumLib
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

    virtual void setConvergenceCheck(IConvergenceCheck &checker) = 0;

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(const std::vector<ICoupledSystem*> &subproblems, UnnamedParameterSet &vars, const ParameterProblemMappingTable &mapping) = 0;
};

}
