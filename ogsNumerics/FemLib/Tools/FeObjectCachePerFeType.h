/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectCachePerFeType.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <map>

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/FemElementFactory.h"


namespace FemLib
{

/**
 * \brief Cache system for finite element objects 
 */
class FeObjectCachePerFeType
{
public:
    explicit FeObjectCachePerFeType(MeshLib::IMesh &msh) : _msh(&msh) {};

    virtual ~FeObjectCachePerFeType()
    {
        BaseLib::releaseObjectsInStdMap(_mapFeObj);
    }

    IFiniteElement* getFeObject(FiniteElementType::type fe_type);

    MeshLib::IMesh* getMesh() const {return _msh;};

private:
    std::map<FiniteElementType::type, IFiniteElement*> _mapFeObj;
    MeshLib::IMesh *_msh;

    DISALLOW_COPY_AND_ASSIGN(FeObjectCachePerFeType);
};

}
