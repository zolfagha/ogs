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
#include <set>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"
#include "MeshLib/Core/IMesh.h"

#include "IDDCGlobalLocalMapping.h"
#include "DDCSlaveDomain.h"

namespace DiscreteLib
{

class DDCSubDomain
{
public:
    DDCSubDomain(MeshLib::IMesh &msh, IDDCGlobaLocalMapping &mapping, std::set<size_t>* list_ghosts=0) : _local_msh(&msh), _map_global_local(&mapping)
    {
        _dom_id = 0;
        if (list_ghosts) {
            _list_ghosts.insert(list_ghosts->begin(), list_ghosts->end());
        }
    }
    virtual ~DDCSubDomain()
    {
        BaseLib::releaseObjectsInStdVector(_list_slaves);
    }

    void setDomainID(size_t i) {_dom_id = i;};
    size_t getDomainID() const {return _dom_id;};

    MeshLib::IMesh* getLoalMesh() {return _local_msh;};

    IDDCGlobaLocalMapping* getGlobalLocalIdMap() {return _map_global_local;};

    size_t getNumberOfGhosts() const {return _list_ghosts.size();};
    std::set<size_t>* getGhostList() {return &_list_ghosts;};

    size_t getNumberOfSlaves() const {return _list_slaves.size();};
    DDCSlaveDomain* getSlaveDomain(size_t i) {return _list_slaves[i];};
private:
    size_t _dom_id;
    MeshLib::IMesh* _local_msh;
    IDDCGlobaLocalMapping* _map_global_local;
    std::vector<DDCSlaveDomain*> _list_slaves;
    std::set<size_t> _list_ghosts;
};

} //end

