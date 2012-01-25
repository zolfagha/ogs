
#pragma once

#include "ShapeFunction.h"

namespace FemLib
{

class IFemTestFunction
{
public:
    virtual MathLib::Matrix<double>* computeTestFunction(MathLib::Matrix<double>& base) = 0;
    virtual MathLib::Matrix<double>* computeGradTestFunction(MathLib::Matrix<double>& dBase) = 0;
};

class FemTestFunctionBubnov : public IFemTestFunction
{
private:
    IFemShapeFunction *_base;

public:
    MathLib::Matrix<double>* computeTestFunction(MathLib::Matrix<double>& base)
    {
        return &base;
    }
    
    MathLib::Matrix<double>* computeGradTestFunction(MathLib::Matrix<double>& dBase)
    {
        return &dBase;
    }

};

}
