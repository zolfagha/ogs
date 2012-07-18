
#pragma once

#include "GeoLib/Core/Point.h"

namespace MeshLib
{
/**
 * \brief 
 *
 * 
 */  
class INode
{
public:
    INode () {};
    virtual ~INode (){};

    virtual size_t getNodeID() const = 0;
    virtual void setNodeID(const size_t &id) = 0;
    virtual const GeoLib::Point* getData() const = 0;
    virtual void setX(const GeoLib::Point &x) = 0; 
};

} // end namespace

