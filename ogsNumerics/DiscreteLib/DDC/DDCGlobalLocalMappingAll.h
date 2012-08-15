/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DomainDecomposition.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/BidirectionalMap.h"
#include "IDDCGlobalLocalMapping.h"

namespace DiscreteLib
{

class DDCGlobaLocalMappingAll : public IDDCGlobaLocalMapping
{
public:
    DDCGlobaLocalMappingAll(BaseLib::BidirectionalMap<size_t, size_t> &map) : _map_global_local(&map)
    {
    }
    virtual ~DDCGlobaLocalMappingAll()
    {
        BaseLib::releaseObject(_map_global_local);
    }
    bool hasGlobal(size_t global_id) 
    {
        return _map_global_local->countInA(global_id)>0;
    }
    size_t global2local(size_t global)
    {
        return _map_global_local->mapAtoB(global);
    }
    size_t local2global(size_t local)
    {
        return _map_global_local->mapBtoA(local);
    }
private:
    BaseLib::BidirectionalMap<size_t, size_t>* _map_global_local;
};

} //end

