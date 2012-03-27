
#pragma once

#include <vector>

#include "MathLib/Function/IFunction.h"

namespace NumLib
{
class ICoupledSystem;
class VariableContainer;
class VariableMappingTable;

/**
 * \brief Interface class for partitioned algorithm
 */
class IPartitionedAlgorithm
{
public:
	typedef MathLib::IFunction Variable;
	typedef VariableContainer MyNamedVariableContainer;
	typedef ICoupledSystem MyCoupledSystem;
	typedef VariableMappingTable MyVariableMappingTable;

	virtual ~IPartitionedAlgorithm() {};

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<MyCoupledSystem*> &subproblems, MyNamedVariableContainer &vars, MyVariableMappingTable &mapping) = 0;
};

}
