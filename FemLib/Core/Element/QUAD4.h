#pragma once

#include "FemLib/Core/IFemElement.h"

namespace FemLib
{

typedef FeBaseIsoparametric<FiniteElementType::QUAD4, 4, FemShapeQuad4, FemIntegrationGaussQuad> QUAD4;

//class QUAD4 : public FeBaseIsoparametric<FiniteElementType::QUAD4, 4, FemShapeQuad4, FemIntegrationGaussQuad>
//{
//protected:
//    virtual IFemShapeFunction* createShapeFunction() const 
//    {
//        return new FemShapeQuad4();
//    }
//
//    virtual IFemIntegration* createIntegrationMethod() const
//    {
//        return new FemIntegrationGaussQuad();
//    }
//};

}
