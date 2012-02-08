
#include "Topology.h"

#include <algorithm>

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace MeshLib
{

TopologySequentialNodes2Elements::TopologySequentialNodes2Elements(const IMesh &msh) 
{
    _node2conn_eles.resize(msh.getNumberOfNodes());
    const size_t nr_ele = msh.getNumberOfElements();
    for (size_t i=0; i<nr_ele; i++) {
        const IElement* e = msh.getElemenet(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            _node2conn_eles[e->getNodeID(j)].push_back(i);
        }
    }
}

TopologyRandomNodes2Elements2::TopologyRandomNodes2Elements2(const IMesh &msh, const std::vector<size_t> &nodes) 
{
    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        IElement* e = msh.getElemenet(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (std::find(nodes.begin(), nodes.end(), e->getNodeID(j))!=nodes.end()) {
                _map_node2conn_eles[e->getID()].push_back(i);
                break;
            }
        }
    }
}

TopologyNode2Nodes::TopologyNode2Nodes(IMesh *msh) 
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


} // end namespace

