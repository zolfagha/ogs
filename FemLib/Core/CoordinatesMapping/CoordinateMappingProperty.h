
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"

namespace FemLib
{

/**
 * \brief Properties of coordinate mapping at particular locations 
 */
struct CoordinateMappingProperty
{
	typedef MathLib::Matrix<double> MatrixType;
    ///// mapping locations in computed coordinates
    //double r[3];
    /// shape function N(r)
	MatrixType *shape_r;
    /// gradient of shape functions, dN(r)/dr
	MatrixType *dshape_dr;
    /// gradient of shape functions, dN(r)/dx
	MatrixType *dshape_dx;
    /// Jacobian matrix, J=dx/dr
	MatrixType *jacobian_dxdr;
    /// determinant of the Jacobian
    double det_jacobian;
    /// inverse of the Jacobian
    MatrixType *inv_jacobian_drdx;

    CoordinateMappingProperty()
    {
        shape_r = new MatrixType();
        dshape_dr = new MatrixType();
        dshape_dx = new MatrixType();
        jacobian_dxdr = new MatrixType();
        inv_jacobian_drdx = new MatrixType();
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
