/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file INode.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "GeoLib/Point.h"

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

