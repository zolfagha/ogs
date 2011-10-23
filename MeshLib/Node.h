
#pragma once

//#include <cmath>
//#include <vector>
#include "Point.h"

//------------------------------------------------------------------------
namespace MeshLib
{

//-----------------------------------------------------------------------------
// Node
//-----------------------------------------------------------------------------
class Node : public GEOLIB::Point
{
private:
    size_t _node_id;
public:
    Node (size_t id, double x, double y, double z) {
        this->_node_id = id;
        this->_x[0] = x;
        this->_x[1] = y;
        this->_x[2] = z;
    };
};

} // end namespace

