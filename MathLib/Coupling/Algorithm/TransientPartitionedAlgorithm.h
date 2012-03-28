
#pragma once

#include <vector>

namespace MathLib
{
class ICoupledSystem;
class ParameterTable;
struct ParameterProblemMappingTable;

/**
 * \brief
 */
class ITransientPartitionedAlgorithm
{
public:
	typedef MathLib::ParameterTable MyNamedVariableContainer;
	typedef MathLib::IFunction Variable;
	typedef MathLib::ICoupledSystem MyCoupledSystem;
	typedef MathLib::ParameterProblemMappingTable MyVariableMappingTable;

	virtual ~ITransientPartitionedAlgorithm() {};

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<MyCoupledSystem*> &subproblems, MyNamedVariableContainer &vars_t_n, MyNamedVariableContainer &vars_t_n1, MyVariableMappingTable &mapping) = 0;
};


}
