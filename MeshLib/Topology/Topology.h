
#pragma once

#include <vector>
#include <set>

#include "MeshLib/Core/IMesh.h"

namespace MeshLib
{

/**
 * \brief Mesh topological data: node <-> connected elements
 */    
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
    std::vector<std::vector<size_t> > _node2conn_eles;
};

/**
 * \brief Mesh topological data: node <-> connected nodes
 */
class TopologyNode2Nodes
{
public:
    TopologyNode2Nodes(IMesh *msh) 
    {
        _node2conn_nodes.resize(msh->getNumberOfNodes());
        const size_t nr_ele = msh->getNumberOfElements();
        for (size_t i=0; i<nr_ele; i++) {
            const IElement* e = msh->getElemenet(i);
            const size_t e_nnodes = e->getNumberOfNodes();
            for (size_t j=0; j<e_nnodes; j++) {
                std::set<size_t> &set_conn_nodes = _node2conn_nodes[e->getNodeID(j)];
                size_t id1 = (j-1+e_nnodes)%e_nnodes;
                size_t id2 = (j+1)%e_nnodes;
                set_conn_nodes.insert(e->getNodeID(id1));
                set_conn_nodes.insert(e->getNodeID(id2));
            }
        }
    }

    const std::set<size_t>& getConnectedNodes(size_t node_id) const 
    {
        return _node2conn_nodes[node_id];
    }

    size_t getNumberOfNodes() const 
    {
        return _node2conn_nodes.size();
    }

private:
    std::vector<std::set<size_t>> _node2conn_nodes;
};

}
