
#pragma once

#include <vector>


namespace NumLib
{
class ICoupledSystem;
class NamedVariableContainer;
class VariableMappingTable;

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
    virtual int solve(std::vector<ICoupledSystem*> &subproblems, NamedVariableContainer &vars, VariableMappingTable &mapping) = 0;
};

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
    virtual int solve(std::vector<ICoupledSystem*> &subproblems, NamedVariableContainer &vars_t_n, NamedVariableContainer &vars_t_n1, VariableMappingTable &mapping) = 0;
};


}
