
#pragma once

#include <vector>

#include "MeshLib/Core/IMesh.h"

namespace MeshLib
{

class TopologyNode2Elements
{
public:
    TopologyNode2Elements(IMesh *msh) {
        _node2conn_eles.resize(msh->getNumberOfNodes());
        const size_t nr_ele = msh->getNumberOfElements();
        for (size_t i=0; i<nr_ele; i++) {
            const IElement* e = msh->getElemenet(i);
            for (size_t j=0; j<e->getNumberOfNodes(); j++) {
                _node2conn_eles[e->getNodeID(j)].push_back(i);
            }
        }
    }

    const std::vector<size_t>& getConnectedElements(size_t node_id) {
        return _node2conn_eles[node_id];
    }
private:
    std::vector<size_t, std::vector<size_t>> _node2conn_eles;
};

}
