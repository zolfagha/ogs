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

#include "IDDCGlobalLocalMapping.h"

namespace DiscreteLib
{

class DDCGlobaLocalMappingOffset : public IDDCGlobaLocalMapping
{
public:
    DDCGlobaLocalMappingOffset(size_t i_start, size_t i_end, size_t offset) : _offset(offset)
    {
        _i_start = i_start;
        _i_end = i_end;
    }
    virtual ~DDCGlobaLocalMappingOffset()
    {
    }
    bool hasGlobal(size_t global_id) 
    {
        return _i_start<=global_id && global_id<_i_end;
    }
    size_t global2local(size_t global)
    {
        return global - _offset;
    }
    size_t local2global(size_t local)
    {
        return local + _offset;
    }
private:
    size_t _i_start;
    size_t _i_end;
    size_t _offset;
};

} //end

