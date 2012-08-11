/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Topology.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
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

class TopologyRandomNodes2Elements
{
public:
    TopologyRandomNodes2Elements(const IMesh &msh, const std::vector<size_t> &nodes);

    const std::vector<size_t>* getConnectedElements(size_t node_id) const {
        std::map<size_t, std::vector<size_t> >::const_iterator itr = _map_node2conn_eles.find(node_id);
        if (itr!=_map_node2conn_eles.end())
            return &itr->second;
        return 0;
    }

private:
    std::map<size_t, std::vector<size_t> > _map_node2conn_eles;
};

class ITopologyNode2Nodes
{
public:
    virtual ~ITopologyNode2Nodes() {};

    virtual const std::set<size_t>& getConnectedNodes(size_t node_id) const = 0;
    virtual size_t getNumberOfNodes() const = 0;
};

/**
 * \brief Mesh topological data: node <-> connected nodes
 */
class TopologyNode2NodesConnectedByEdges : public ITopologyNode2Nodes
{
public:
    TopologyNode2NodesConnectedByEdges(IMesh *msh);
    virtual ~TopologyNode2NodesConnectedByEdges() {};

    const std::set<size_t>& getConnectedNodes(size_t node_id) const 
    {
        return _node2conn_nodes[node_id];
    }

    size_t getNumberOfNodes() const 
    {
        return _node2conn_nodes.size();
    }

private:
    std::vector<std::set<size_t> > _node2conn_nodes;
};

class TopologyNode2NodesConnectedByElements : public ITopologyNode2Nodes
{
public:
    TopologyNode2NodesConnectedByElements(IMesh *msh);
    virtual ~TopologyNode2NodesConnectedByElements() {};

    const std::set<size_t>& getConnectedNodes(size_t node_id) const 
    {
        return _node2conn_nodes[node_id];
    }

    size_t getNumberOfNodes() const 
    {
        return _node2conn_nodes.size();
    }

private:
    std::vector<std::set<size_t> > _node2conn_nodes;
};


}
