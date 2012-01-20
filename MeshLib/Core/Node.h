
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
template<typename Tpos>
class Node : public INode<Tpos>
{
private:
    size_t _node_id;
    std::set<size_t> _connected_nodes;
    Tpos _x;
public:
    Node (size_t id, const Tpos &x) {
        this->_node_id = id;
        this->_x = x;
    };

    size_t getNodeID(size_t id) { return _node_id;};
    void setNodeID(size_t id) { _node_id = id;};
    virtual const Tpos* getData() const 
    {
        return &_x;
    };
    void setX(const Tpos &x) 
    {
        _x = x;
    };
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

