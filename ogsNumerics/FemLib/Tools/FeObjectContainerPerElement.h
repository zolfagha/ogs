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

#include <map>

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Core/Element/FemElementFactory.h"
#include "IFeObjectContainer.h"

namespace FemLib
{

/**
 * \brief Finite element object containers
 */
class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    FeObjectContainerPerElement(MeshLib::IMesh &msh, size_t ele_size) : _msh(&msh)
    {
        _vec_fem.resize(ele_size);
    }
    virtual ~FeObjectContainerPerElement()
    {
        BaseLib::releaseObjectsInStdVector(_vec_fem);
    }

    void addFiniteElement(size_t i, FiniteElementType::type fe_type)
    {
        _vec_fem[i] = FemElementFactory::create(fe_type, *_msh);
    }

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e) 
    {
        return _vec_fem[e.getID()];
    }

private:
    MeshLib::IMesh* _msh;
    std::vector<IFiniteElement*> _vec_fem;
    DISALLOW_COPY_AND_ASSIGN(FeObjectContainerPerElement);
};

}
