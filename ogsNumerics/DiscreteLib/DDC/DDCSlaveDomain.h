/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DDCSlaveDomain.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

namespace DiscreteLib
{

class DDCSlaveDomain
{
public:
    size_t getDomainID() const {return _dom_id;};
    size_t getNumberOfSlaveObjects() const {return _list_slave_objects.size();};
    size_t getSlaveObjectID(size_t i) const {return _list_slave_objects[i];};
private:
    size_t _dom_id;
    std::vector<size_t> _list_slave_objects;
};

} //end

