
#pragma once

#include <limits>

#include "GeoLib/AxisAlignedBoundingBox.h"
#include "CoordinateSystem.h"

namespace MeshLib
{

class MeshGeometricProperty
{
public:
    MeshGeometricProperty()
    {
        _min_edge_len = std::numeric_limits<double>::epsilon();
    }
    double getMinEdgeLength() const
    {
        return _min_edge_len;
    }

    void setMinEdgeLength(double l)
    {
        _min_edge_len = l;
    }

    const CoordinateSystem* getCoordinateSystem() const
    {
        return &_coord;
    }
    void setCoordinateSystem(const CoordinateSystemType::type &coord)
    {
        _coord.setType(coord);
    }

    const GeoLib::AxisAlignedBoundingBox& getBoundingBox() const
    {
        return _aabox;
    }

    GeoLib::AxisAlignedBoundingBox& getBoundingBox()
    {
        return _aabox;
    }

private:
    double _min_edge_len;
    CoordinateSystem _coord;
    GeoLib::AxisAlignedBoundingBox _aabox;
};

}// end namespace

