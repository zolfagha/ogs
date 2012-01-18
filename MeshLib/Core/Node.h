
#pragma once

//#include <cmath>
#include <set>
#include "INode.h"
#include "GeoLib/Point.h"

//------------------------------------------------------------------------
namespace MeshLib
{

//-----------------------------------------------------------------------------
// Node
//-----------------------------------------------------------------------------
class Node : public INode, public GeoLib::Point
{
private:
    size_t _node_id;
    std::set<size_t> _connected_nodes;
public:
    Node (size_t id, double x, double y, double z) {
        this->_node_id = id;
        this->_x[0] = x;
        this->_x[1] = y;
        this->_x[2] = z;
    };

    size_t getNodeID(size_t id) { return _node_id;};
    void setNodeID(size_t id) { _node_id = id;};

    const std::set<size_t>& getConnectedNodes() const
    {
        return _connected_nodes;
    }

    void addConnectedNode(size_t node_id)
    {
        _connected_nodes.insert(node_id);
    }

};

} // end namespace

