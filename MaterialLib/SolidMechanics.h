
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"

namespace SolidMechanics
{

class TransformatiomMatrix
{
    void set2D(const MathLib::Matrix<double> &_matVec)
    {
        MathLib::Matrix<double> _matTensor(4,4);
        _matTensor(0,0) = _matVec(0,0)*_matVec(0,0);
        _matTensor(0,1) = _matVec(1,0)*_matVec(1,0);
        _matTensor(0,2) = 0.0;
        _matTensor(0,3) = -_matVec(0,0)*_matVec(1,0);
        _matTensor(1,0) = _matVec(0,1)*_matVec(0,1);
        _matTensor(1,1) = _matVec(1,1)*_matVec(1,1);
        _matTensor(1,2) = 0.0;
        _matTensor(1,3) = -_matVec(0,1)*_matVec(1,1);
        _matTensor(2,0) = 0.0;
        _matTensor(2,1) = 0.0;
        _matTensor(2,2) = 1.0;
        _matTensor(2,3) = 0.0;
        _matTensor(3,0) = -2.0*_matVec(0,0)*_matVec(0,1);
        _matTensor(3,1) = -2.0*_matVec(1,0)*_matVec(1,1);
        _matTensor(3,2) = 0.0;
        _matTensor(3,3) = _matVec(0,0)*_matVec(1,1)+_matVec(1,0)*_matVec(0,1);
    };

    void set3D()
    {

    };
};

}
