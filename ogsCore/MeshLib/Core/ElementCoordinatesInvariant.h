/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementCoordinatesInvariant
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "GeoLib/Point.h"

#include "MeshLib/Core/IElementCoordinatesMapping.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/CoordinateSystem.h"

namespace MeshLib
{

/**
 * \brief EleMapInvariant keep original coordinates.
 */
class ElementCoordinatesInvariant : public IElementCoordinatesMapping
{
public:
    ///
    ElementCoordinatesInvariant(const IMesh* msh, IElement* e) : _msh(msh), _e(e) 
    {
    };
    
    ///
    virtual ~ElementCoordinatesInvariant() {};

    ///
    GeoLib::Point* getNodePoint(size_t local_id)
    {
        return (GeoLib::Point*)_msh->getNodeCoordinatesRef(_e->getNodeID(local_id));
    }

private:
    const IMesh* _msh;
    IElement* _e;
};

}
