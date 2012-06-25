
#pragma once

#include <iostream>

#include "MeshLib/Core/IMesh.h"
#include "CoordinateMappingProperty.h"

namespace FemLib
{

/**
 * \brief IFemCoordinatesMapping is an interface class for geometrical mapping between actual coordinates and computational coordinates.
 */
class IFemCoordinatesMapping
{
public:
	virtual ~IFemCoordinatesMapping() {};
    /// initialize element
    virtual void initialize(MeshLib::IElement &ele) = 0;
    /// compute shape functions
    virtual const CoordinateMappingProperty* compute(const double* natural_pt) = 0;
    /// map natural coordinates to physical coordinates
    virtual void mapToPhysicalCoordinates(const double* natural_pt, double*) = 0;
    /// map physical coordinates to natural coordinates
    virtual void mapFromPhysicalCoordinates(const double* physical_pt, double*) = 0;

};

}
