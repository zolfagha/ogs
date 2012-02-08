
#pragma once

#include <vector>
#include <map>
#include <set>

namespace MeshLib
{
class IMesh;

/**
 * \brief Mesh topological data: node <-> connected elements
 */    
class TopologySequentialNodes2Elements
{
public:
    TopologySequentialNodes2Elements(const IMesh &msh);

    const std::vector<size_t>& getConnectedElements(size_t node_id) const {
        return _node2conn_eles[node_id];
    }

private:
    std::vector<std::vector<size_t> > _node2conn_eles;
};

class TopologyRandomNodes2Elements2
{
public:
    TopologyRandomNodes2Elements2(const IMesh &msh, const std::vector<size_t> &nodes);

    const std::vector<size_t>* getConnectedElements(size_t node_id) const {
        std::map<size_t, std::vector<size_t> >::const_iterator itr = _map_node2conn_eles.find(node_id);
        if (itr!=_map_node2conn_eles.end())
            return &itr->second;
        return 0;
    }

private:
    std::map<size_t, std::vector<size_t> > _map_node2conn_eles;
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
