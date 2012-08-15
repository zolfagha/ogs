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

#include "BaseLib/CodingTools.h"

#include "DdcEnums.h"
#include "DDCSubdDomain.h"

namespace DiscreteLib
{

class DDCGlobal
{
public:
    DDCGlobal(DecompositionType::type ddc_type) : _ddc_type(ddc_type), _id(0)
    {
        _n_discrete_pt = 0;
    };
    virtual ~DDCGlobal()
    {
        BaseLib::releaseObjectsInStdVector(_list_dom);
    }

    size_t getID() const {return _id;};
    void setID(size_t i) {_id = i;};

    DecompositionType::type getDecompositionType() const {return _ddc_type;};

    size_t addSubDomain(DDCSubDomain* sub)
    {
        _list_dom.push_back(sub);
        sub->setDomainID(_list_dom.size()-1);
        if (_ddc_type==DecompositionType::Node) {
            _n_discrete_pt += sub->getLoalMesh()->getNumberOfNodes();
        } else {
            _n_discrete_pt += sub->getLoalMesh()->getNumberOfElements();
        }
        _n_discrete_pt -= sub->getNumberOfGhosts();
        return sub->getDomainID();
    }

    size_t getNumberOfSubDomains() const {return _list_dom.size();};
    DDCSubDomain* getSubDomain(size_t i) {return _list_dom[i];};

    size_t findDomainID(size_t global_obj_id) const
    {
        for (size_t i=0; i<_list_dom_start_obj.size(); ++i) {
            if (_list_dom_start_obj[i] <= global_obj_id)
                return i;
        }
        return 0; //error
    }

    size_t getTotalNumberOfDecomposedObjects() const {return _n_discrete_pt;}; 
private:
    DISALLOW_COPY_AND_ASSIGN(DDCGlobal);

private:
    DecompositionType::type _ddc_type;
    std::vector<size_t> _list_dom_start_obj;
    std::vector<DDCSubDomain*> _list_dom;
    size_t _n_discrete_pt;
    size_t _id;
};

} //end

