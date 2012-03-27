
#pragma once

#include <vector>


namespace NumLib
{
class ICoupledSystem;
class VariableContainer;
class VariableMappingTable;

/**
 * \brief
 */
class ITransientPartitionedAlgorithm
{
public:
	typedef VariableContainer MyNamedVariableContainer;
	typedef MathLib::IFunction Variable;
	typedef ICoupledSystem MyCoupledSystem;
	typedef VariableMappingTable MyVariableMappingTable;

	virtual ~ITransientPartitionedAlgorithm() {};
    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<MyCoupledSystem*> &subproblems, MyNamedVariableContainer &vars_t_n, MyNamedVariableContainer &vars_t_n1, MyVariableMappingTable &mapping) = 0;
};


}
