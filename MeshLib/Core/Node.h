
#pragma once

//#include <cmath>
#include <set>

#include "MathLib/Vector.h"
#include "GeoLib/Core/Point.h"

#include "INode.h"

namespace MeshLib
{

/**
 * 
 */
class Node : public INode
{
private:
    size_t _node_id;
    std::set<size_t> _connected_nodes;
    GeoLib::Point _x;
public:
    Node() {};
    Node (size_t id, const GeoLib::Point &x) {
        this->_node_id = id;
        this->_x = x;
    };

    size_t getNodeID() const { return _node_id;};
    void setNodeID(const size_t &id) { _node_id = id;};
    virtual const GeoLib::Point* getData() const 
    {
        return &_x;
    };
    void setX(const GeoLib::Point &x) 
    {
        _x = x;
    };
};

} // end namespace

