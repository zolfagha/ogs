
#pragma once

#include <vector>

namespace NumLib
{

class ICoupledSystem;

/**
 * \brief Mapping between variables
 */
struct VariableMappingTable
{
	typedef MathLib::IFunction Variable;
    typedef std::pair<ICoupledSystem*,size_t> PairSysVarId;
    typedef std::pair<size_t, size_t> PairInputVar;
    typedef std::vector<PairInputVar> ListOfInputVar;
    /// List of input variables
    std::vector<ListOfInputVar> _list_subproblem_input_source;
    /// Parameter Id
    std::vector<PairSysVarId> _map_paraId2subproblem;
};


}
