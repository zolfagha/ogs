
#pragma once

namespace MeshLib
{

struct CoordinateSystemType
{
    enum type {
        X = 10,
        Y = 11,
        Z = 12,
        XY = 21,
        XZ = 22,
        YZ = 23,
        XYZ = 32,
        INVALID = -1
    };
};

/**
 * \brief Coordinate systems
 *
 *
 */
class CoordinateSystem
{
public:

    CoordinateSystem() {
        _type = CoordinateSystemType::INVALID;
    }

    CoordinateSystem(CoordinateSystemType::type coord) {
        _type = coord;
    }

    void setType(CoordinateSystemType::type coord) {
        _type = coord;
    }

    /// get this coordinate type
    CoordinateSystemType::type getType() const {
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
    CoordinateSystemType::type _type; 
};

}
