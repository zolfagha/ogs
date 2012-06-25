
#pragma once

#include "FemLib/Core/DataType.h"

namespace FemLib
{

/**
 * \brief Properties of coordinate mapping at particular locations 
 */
struct CoordinateMappingProperty
{
    ///// mapping locations in computed coordinates
    //double r[3];
    /// shape function N(r)
	LocalMatrix *shape_r;
    /// gradient of shape functions, dN(r)/dr
	LocalMatrix *dshape_dr;
    /// gradient of shape functions, dN(r)/dx
	LocalMatrix *dshape_dx;
    /// Jacobian matrix, J=dx/dr
	LocalMatrix *jacobian_dxdr;
    /// determinant of the Jacobian
    double det_jacobian;
    /// inverse of the Jacobian
    LocalMatrix *inv_jacobian_drdx;

    CoordinateMappingProperty()
    {
        shape_r = new LocalMatrix();
        dshape_dr = new LocalMatrix();
        dshape_dx = new LocalMatrix();
        jacobian_dxdr = new LocalMatrix();
        inv_jacobian_drdx = new LocalMatrix();
    }

    ~CoordinateMappingProperty()
    {
        delete shape_r;
        delete dshape_dr;
        delete dshape_dx;
        delete jacobian_dxdr;
        delete inv_jacobian_drdx;
    }
};

}
