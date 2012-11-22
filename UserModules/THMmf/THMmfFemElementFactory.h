/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file THMmfFemElementFactory.h
 *
 * Created on 2012-11-22 by Norihiro Watanabe
 */

#pragma once

#include "THMmfFiniteElementType.h"

namespace MeshLib
{
class IMesh;
}
namespace FemLib
{
class IFiniteElement;
}


namespace THMmf
{

class THMmfFemElementFactory
{
public:
    static FemLib::IFiniteElement* create(THMmfFiniteElementType::type fe_type, MeshLib::IMesh* msh);
};

}
