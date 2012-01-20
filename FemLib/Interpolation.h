
#pragma once

#include "ShapeFunction.h"

namespace FemLib
{

class IFemInterpolation
{
public:
    virtual void computeX() = 0;
    virtual void computeDX() = 0;
private:
    IFemShapeFunction *_shape;
};


}
