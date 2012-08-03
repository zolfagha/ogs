/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectCachePerFeType.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FeObjectCachePerFeType.h"

namespace FemLib
{

IFiniteElement* FeObjectCachePerFeType::getFeObject(FiniteElementType::type fe_type)
{
    IFiniteElement *fe = 0;
    if (_mapFeObj.count(fe_type)==0) {
        fe = FemElementFactory::create(fe_type, *_msh);
        _mapFeObj[fe_type] = fe;
    } else {
        fe = _mapFeObj[fe_type];
    }
    return fe;
}

}
