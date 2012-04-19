
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Mapping.h"
#include "FemLib/Core/ShapeFunction/IFemShapeFunction.h"
#include "IFemCoordinatesMapping.h"

namespace FemLib
{

/**
 * \brief Mapping element shapes to natural coordinates
 *
 * FemNaturalCoordinates mapping converts element shapes in physical coordinates (x,y,z) to that in natural coordinates (r,s,t).
 * - Given physical coordinates should correspond to dimensions of the element, i.e (x,y) for triangles
 * - x(r,s,t) = N(r,s,t) x_i
 */
class FemNaturalCoordinates : public IFemCoordinatesMapping
{
public:
    FemNaturalCoordinates(IFemShapeFunction *shape) 
    {
        _shape = shape;
        _prop = new CoordinateMappingProperty();
    }

    virtual ~FemNaturalCoordinates()
    {
        delete _prop;
    }

    /// initialize element
    virtual void initialize(MeshLib::IElement &ele);

    virtual const CoordinateMappingProperty* getProperties() const {return _prop;};

    /// compute mapping properties at the given location in natural coordinates
    virtual const CoordinateMappingProperty* compute(const double* natural_pt);

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param natural_pt
    /// @return physical_pt
    void mapToPhysicalCoordinates(const double* natural_pt, double* physical_pt);

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param natural_pt
    /// @return physical_pt
    void mapToPhysicalCoordinates(const CoordinateMappingProperty* prop, double* physical_pt);

    /// compute natural coordinates at the given natural coordinates.
    /// Assuming \f$ r=0 \f$ at \f$ x = \bar{x}_{avg} \f$, natural coordinates can be calculated as
    /// \f[
    ///    \mathbf{r} = (\mathbf{J}^-1)^T * (\mathbf{x} - \bar{\mathbf{x}}_{avg})
    /// \f]
    ///
    /// @param physical_pt
    /// @return natural_pt
    void mapFromPhysicalCoordinates(const double* physical_pt, double* natural_pt);

private:
    CoordinateMappingProperty* _prop;
    IFemShapeFunction* _shape;
    MeshLib::IElement* _ele;
    MathLib::Matrix<double> x;

};


}
