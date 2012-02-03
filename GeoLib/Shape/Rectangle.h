
#pragma once

#include <vector>

#include "GeoLib/Core/Point.h"
#include "GeoLib/Core/Polyline.h"
#include "Base/MemoryTools.h"

namespace GeoLib
{

class Rectangle
{
public:
    Rectangle(const Point &pt_lower_left, const Point &pt_upper_right) {
        const double z = pt_lower_left[2];
        _pnt_vec.push_back(new Point(pt_lower_left));
        _pnt_vec.push_back(new Point(pt_upper_right[0], pt_lower_left[1], z));
        _pnt_vec.push_back(new Point(pt_upper_right));
        _pnt_vec.push_back(new Point(pt_lower_left[0], pt_upper_right[1], z));
    };
    virtual ~Rectangle() {
        Base::destroyStdVectorWithPointers(_pnt_vec);
        Base::destroyStdVectorWithPointers(_poly_vec);
    }

    Polyline* getDomainPolyline() 
    {
        Polyline *poly = new Polyline(_pnt_vec);
        poly->addPoint(0);
        poly->addPoint(1);
        poly->addPoint(2);
        poly->addPoint(3);
        _poly_vec.push_back(poly);
        return poly;
    };

    Polyline* getLeft()
    {
        Polyline *poly = new Polyline(_pnt_vec);
        poly->addPoint(0);
        poly->addPoint(3);
        _poly_vec.push_back(poly);
        return poly;
    };

    Polyline* getRight() 
    {
        Polyline *poly = new Polyline(_pnt_vec);
        poly->addPoint(1);
        poly->addPoint(2);
        _poly_vec.push_back(poly);
        return poly;
    };

    Polyline* getTop() 
    {
        Polyline *poly = new Polyline(_pnt_vec);
        poly->addPoint(2);
        poly->addPoint(3);
        _poly_vec.push_back(poly);
        return poly;
    };

    Polyline* getBottom() 
    {
        Polyline *poly = new Polyline(_pnt_vec);
        poly->addPoint(0);
        poly->addPoint(1);
        _poly_vec.push_back(poly);
        return poly;
    };

private:
    std::vector<Point*> _pnt_vec;
    std::vector<Polyline*> _poly_vec;
};

}
