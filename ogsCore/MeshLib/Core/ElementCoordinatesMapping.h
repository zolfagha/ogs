/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementCoordinatesMapping.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "GeoLib/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/CoordinateSystem.h"

namespace MeshLib
{

/**
 * \brief Interface of methods mapping coordinates of elements.
 *
 *
 */
class IElementCoordinatesMapping
{
public:
    virtual ~IElementCoordinatesMapping() {};
    /// get mapped coordinates of nodes
    virtual GeoLib::Point* getNodePoint(size_t node_id) = 0;
};

/**
 * \brief EleMapInvariant keep original coordinates.
 */
class EleMapInvariant : public IElementCoordinatesMapping
{
public:
    EleMapInvariant(const IMesh &msh, IElement &e) 
    {
        _msh = &msh;
        _e = &e;
    };
    virtual ~EleMapInvariant() {};

    GeoLib::Point* getNodePoint(size_t local_id) 
    {
        return (GeoLib::Point*)_msh->getNodeCoordinatesRef(_e->getNodeID(local_id));
    }

private:
    const IMesh* _msh;
    IElement *_e;
};

/**
 * \brief Mapping local coordinates of elements.
 *
 * 
 */
class EleMapLocalCoordinates : public IElementCoordinatesMapping 
{
public:
    ///
    EleMapLocalCoordinates(const IMesh &msh, IElement &e, const CoordinateSystem &coordinate_system) : _matR2original(0)
    {
        assert (e.getDimension() <= coordinate_system.getDimension());
        _msh = &msh;

        //if (e->getDimension()==coordinate_system->getDimension()) {
        //    flip(e, coordinate_system);
        ////} else if (e->getDimension() < coordinate_system->getDimension()) {
        ////    rotate(e, coordinate_system);
        //}
        rotate(e, coordinate_system);
    };

    virtual ~EleMapLocalCoordinates()
    {
        BaseLib::releaseObjectsInStdVector(_point_vec);
        if (_matR2original!=0)
            delete _matR2original;
    }

    virtual GeoLib::Point* getNodePoint(size_t node_id) {
        return _point_vec[node_id];
    }

private:
    const IMesh* _msh;
    std::vector<GeoLib::Point*> _point_vec;
    MathLib::Matrix<double> *_matR2original;
    MathLib::Matrix<double> *_matR2local;

    ///
    void flip(IElement &e, const CoordinateSystem &coordinate_system);
    ///
    void rotate(IElement &e, const CoordinateSystem &coordinate_system);
    // x=Rx' where x is original coordinates and x' is local coordinates
    void getRotationMatrixToOriginal(const IElement &e, const CoordinateSystem &coordinate_system, const std::vector<GeoLib::Point> &vec_pt);
};

}
