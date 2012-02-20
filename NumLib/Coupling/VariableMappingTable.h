
#pragma once

#include <vector>

namespace NumLib
{

class ICoupledSystem;

class VariableMappingTable
{
public:
    typedef std::pair<ICoupledSystem*,size_t> PairSysVarId;
    typedef std::pair<size_t, size_t> PairInputVar;
    typedef std::vector<PairInputVar> ListOfInputVar;
    std::vector<ListOfInputVar> _list_subproblem_input_source;
    std::vector<PairSysVarId> _map_paraId2subproblem;
};


}