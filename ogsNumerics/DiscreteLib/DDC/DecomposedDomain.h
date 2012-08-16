/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file DecomposedDomain.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

#include "DdcEnums.h"

namespace DiscreteLib
{

class SubDomain;

/**
 * \brief Decomposed domain class
 *
 * This class owns
 * - decomposition type
 * - sub-domains
 */
class DecomposedDomain
{
public:
    explicit DecomposedDomain(DecompositionType::type ddc_type)
    : _ddc_type(ddc_type), _n_discrete_pt(0), _id(0)
    { };

    virtual ~DecomposedDomain()
    {
        BaseLib::releaseObjectsInStdVector(_list_dom);
    }

    size_t getID() const {return _id;};

    void setID(size_t i) {_id = i;};

    DecompositionType::type getDecompositionType() const {return _ddc_type;};

    size_t getTotalNumberOfDecomposedObjects() const {return _n_discrete_pt;};

    size_t getNumberOfSubDomains() const {return _list_dom.size();};

    size_t addSubDomain(SubDomain* sub);

    SubDomain* getSubDomain(size_t i) {return _list_dom[i];};

    size_t findDomainID(size_t global_obj_id) const;

private:
    DISALLOW_COPY_AND_ASSIGN(DecomposedDomain);

private:
    DecompositionType::type _ddc_type;
    std::vector<size_t> _list_dom_start_obj;
    std::vector<SubDomain*> _list_dom;
    size_t _n_discrete_pt;
    size_t _id;
};

} //end

