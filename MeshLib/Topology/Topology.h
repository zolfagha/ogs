
#pragma once

#include <vector>
#include <set>

namespace MeshLib
{
class IMesh;

/**
 * \brief Mesh topological data: node <-> connected elements
 */    
class TopologyNode2Elements
{
public:
    TopologyNode2Elements(const IMesh &msh);

    const std::vector<size_t>& getConnectedElements(size_t node_id) const {
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
    TopologyNode2Nodes(IMesh *msh);

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
