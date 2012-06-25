
#pragma once

#include "IFemShapeFunction.h"

namespace FemLib
{

template <size_t N_DIM, size_t N_NODES>
class TemplateShapeFunction : public IFemShapeFunction
{
public:
    virtual void computeShapeFunction(const double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(const double* pt, double* dN) = 0;
    virtual ~TemplateShapeFunction() {};
protected:
};
  
}
