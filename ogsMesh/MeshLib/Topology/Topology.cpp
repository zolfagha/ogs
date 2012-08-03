/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Topology.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

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

TopologyRandomNodes2Elements::TopologyRandomNodes2Elements(const IMesh &msh, const std::vector<size_t> &nodes) 
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

TopologyNode2NodesConnectedByEdges::TopologyNode2NodesConnectedByEdges(IMesh *msh) 
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

TopologyNode2NodesConnectedByElements::TopologyNode2NodesConnectedByElements(IMesh *msh) 
{
    _node2conn_nodes.resize(msh->getNumberOfNodes());
    const size_t nr_ele = msh->getNumberOfElements();
    for (size_t i=0; i<nr_ele; i++) {
        const IElement* e = msh->getElemenet(i);
        const size_t e_nnodes = e->getNumberOfNodes();
        for (size_t j=0; j<e_nnodes; j++) {
            std::set<size_t> &set_conn_nodes = _node2conn_nodes[e->getNodeID(j)];
            for (size_t k=0; k<e_nnodes; k++) {
                if (j==k) continue;
                set_conn_nodes.insert(e->getNodeID(k));
            }
        }
    }
}

} // end namespace

