/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Line.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "GeoLib/Core/Point.h"
#include "GeoLib/Core/Polyline.h"

namespace GeoLib
{

class Line
{
public:
    Line(const Point &pt1, const Point &pt2)
    {
        _pnt_vec.push_back(new Point(pt1));
        _pnt_vec.push_back(new Point(pt2));
    }

    virtual ~Line() {
        BaseLib::releaseObjectsInStdVector(_pnt_vec);
        BaseLib::releaseObjectsInStdVector(_poly_vec);
    }

    Polyline* getPolyline()
    {
        Polyline *poly = new Polyline(_pnt_vec);
        poly->addPoint(0);
        poly->addPoint(1);
        _poly_vec.push_back(poly);
        return poly;
    }

    Point* getPoint1()
    {
        return _pnt_vec[0];
    }

    Point* getPoint2()
    {
        return _pnt_vec[1];
    }

private:
    std::vector<Point*> _pnt_vec;
    std::vector<Polyline*> _poly_vec;
};

}
