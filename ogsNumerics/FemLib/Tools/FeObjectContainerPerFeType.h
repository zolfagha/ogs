/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectContainerPerFeType.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "IFeObjectContainer.h"
#include "FeObjectCachePerFeType.h"
#include "FemElementCatalog.h"
#include "MeshElementShapeToFemElementType.h"

namespace FemLib
{

/**
 * \brief Finite element object containers based on element shape types
 */
class FeObjectContainerPerFeType
 : public IFeObjectContainer
{
public:
    /**
     *
     * @param fe_catalog    FE type catalog
     * @param shape2feType  Link between shape type and FE type
     * @param msh           Mesh
     */
    FeObjectContainerPerFeType(const FemElementCatalog* fe_catalog, const MeshElementShapeToFemElementType* shape2feType, MeshLib::IMesh* msh)
    : _cache(fe_catalog, msh), _shape2feType(shape2feType)
    {
    };

    /**
     * Copy constructor
     * @param src
     */
    FeObjectContainerPerFeType(const FeObjectContainerPerFeType &src)
    : _cache(src._cache), _shape2feType(src._shape2feType)
    {
    }

    /**
     *
     */
    virtual ~FeObjectContainerPerFeType() {};

    /**
     *
     * @param e     Mesh element
     * @return a pointer to IFiniteElement object
     */
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e)
    {
        int fe_type = _shape2feType->getFeType(e, e.getCurrentOrder());
        IFiniteElement* fe = _cache.getFeObject(fe_type);
        fe->configure(*const_cast<MeshLib::IElement*>(&e));
        return fe;
    }

private:
    FeObjectCachePerFeType _cache;
    const MeshElementShapeToFemElementType* _shape2feType;
};



}
