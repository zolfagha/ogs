/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementCoordinatesMappingLocal
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
 * \brief Mapping local coordinates of elements.
 *
 * 
 */
class ElementCoordinatesMappingLocal : public IElementCoordinatesMapping 
{
public:
    ///
    ElementCoordinatesMappingLocal(const IMesh* msh, IElement &e, const CoordinateSystem &coordinate_system)
    : _msh(msh), _matR2original(0)
    {
        assert (e.getDimension() <= coordinate_system.getDimension());

        //if (e->getDimension()==coordinate_system->getDimension()) {
        //    flip(e, coordinate_system);
        ////} else if (e->getDimension() < coordinate_system->getDimension()) {
        ////    rotate(e, coordinate_system);
        //}
        rotate(e, coordinate_system);
    };

    ///
    virtual ~ElementCoordinatesMappingLocal()
    {
        BaseLib::releaseObjectsInStdVector(_point_vec);
        if (_matR2original!=0)
            delete _matR2original;
    }

    ///
    virtual GeoLib::Point* getNodePoint(size_t node_id) 
    {
        return _point_vec[node_id];
    }

private:
    ///
    void flip(IElement &e, const CoordinateSystem &coordinate_system);
    ///
    void rotate(IElement &e, const CoordinateSystem &coordinate_system);
    // x=Rx' where x is original coordinates and x' is local coordinates
    void getRotationMatrixToOriginal(const IElement &e, const CoordinateSystem &coordinate_system, const std::vector<GeoLib::Point> &vec_pt);

private:
    const IMesh* _msh;
    std::vector<GeoLib::Point*> _point_vec;
    MathLib::Matrix<double> *_matR2original;
    MathLib::Matrix<double> *_matR2local;
};

}
