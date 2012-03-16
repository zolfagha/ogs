
#include "DoF.h"


namespace DiscreteLib
{
#if 0
size_t DofMap::setEqsIDSequnetual(size_t eqs_id_begin)
{
    size_t eqs_id = eqs_id_begin;
    for (size_t i=0; i<_map_node_id2eqs_id.size(); i++) {
        if (isActiveDoF(i))
            _map_node_id2eqs_id[i] = eqs_id++;
        else
            _map_node_id2eqs_id[i] = -1;
    }
    return eqs_id;
}
#endif


} //end
