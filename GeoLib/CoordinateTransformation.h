
#pragma once

#include <cmath>
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/MathTools.h"

namespace GeoLib
{

class ICoordinateSystem
{
public:
    virtual size_t getDimension() = 0;
    virtual const MathLib::Matrix<double>* getTransformationMatrix() const = 0;
};

/**
 *  create a rotation matrix with an angle. It rotates points in the xy plane 
 *  counterclockwise through an angle theta about the origin.
 *  \f[
 *      \mathbf{R}_x
 *       = 
 *      \left [ \begin{array}{cc} \cos \theta & -\sin \theta \\ \sin \theta & \cos \theta \end{array} \right ]
 *  \f]
 */
void createLocalCoordinateSystem2D(double rotate_theta, MathLib::Matrix<double> &_matX)
{
    _matX(0,0) = cos(rotate_theta);
    _matX(0,1) = -sin(rotate_theta);
    _matX(1,0) = sin(rotate_theta);
    _matX(1,1) = cos(rotate_theta);
}

/**
 * create a rotation matrix from a normal vector \f$\mathbf{n}\f$
 *  \f[
 *      \mathbf{R}_x
 *       = 
 *      \left [ \begin{array}{cc} n_y & -n_x \\ n_x & n_y \end{array} \right ]
 *  \f]
 */
void createLocalCoordinateSystem2D(double *normal_vec, MathLib::Matrix<double> &_matX) {
    _matX(0,0) = normal_vec[1];
    _matX(0,1) = -normal_vec[0];
    _matX(1,0) = normal_vec[0];
    _matX(1,1) = normal_vec[1];
}


}
