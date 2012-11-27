/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectContainerPerElement.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "IFeObjectContainer.h"
#include "FemElementCatalog.h"

namespace FemLib
{

/**
 * \brief Finite element object containers
 */
class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    /**
     *
     * @param msh
     * @param ele_size
     */
    FeObjectContainerPerElement(const FemElementCatalog* fe_catalog, MeshLib::IMesh* msh, size_t ele_size)
    : _fe_catalog(fe_catalog), _msh(msh)
    {
        _vec_fem.resize(ele_size);
    }

    ///
    virtual ~FeObjectContainerPerElement()
    {
        BaseLib::releaseObjectsInStdVector(_vec_fem);
    }

    /**
     *
     * @param i
     * @param fe_type
     */
    void addFiniteElement(size_t i, int fe_type)
    {
        _vec_fem[i] = _fe_catalog->createFeObject(fe_type, _msh);
    }

    /**
     *
     * @param e
     * @return
     */
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e) 
    {
        return _vec_fem[e.getID()];
    }

private:
    DISALLOW_COPY_AND_ASSIGN(FeObjectContainerPerElement);

private:
    const FemElementCatalog* _fe_catalog;
    MeshLib::IMesh* _msh;
    std::vector<IFiniteElement*> _vec_fem;
};

}
