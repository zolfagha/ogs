
#pragma once

namespace MeshLib
{

/**
 * \brief Coordinate systems
 *
 *
 */
class CoordinateSystem
{
public:
    enum CoordinateSystemType {
        X = 10,
        Y = 11,
        Z = 12,
        XY = 21,
        XZ = 22,
        XYZ = 32
    };

    CoordinateSystem(CoordinateSystemType coord) {
        _type = coord;
    }

    /// get this coordinate type
    CoordinateSystemType getType() const {
        return _type;
    }

    /// get dimension size
    size_t getDimension() const {
        return _type / 10;
    }

    /// has z dimension
    bool hasZ() const {
        return (_type % 10 == 2);
    }

private:
    CoordinateSystemType _type; 
};
}
