
#pragma once


namespace FemLib
{

class IFemShapeFunction
{
public:
    virtual void computeShapeFunction(const double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(const double* pt, double* dN) = 0;
    virtual ~IFemShapeFunction() {};
};

}
