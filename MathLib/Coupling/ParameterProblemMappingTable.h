
#pragma once

#include <vector>

namespace MathLib
{

class ICoupledSystem;

/**
 * \brief Mapping between variables
 */
struct ParameterProblemMappingTable
{
    typedef std::pair<ICoupledSystem*,size_t> PairSysVarId;
    typedef std::pair<size_t, size_t> PairInputVar;
    typedef std::vector<PairInputVar> ListOfInputVar;
    /// List of input variables
    std::vector<ListOfInputVar> _list_subproblem_input_source;
    /// Parameter Id
    std::vector<PairSysVarId> _map_paraId2subproblem;
};


}
