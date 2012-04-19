
#pragma once

#include "FemLib/Core/Extrapolation/FemExtrapolation.h"
#include "FiniteElementType.h"
#include "TemplateIsoparametric.h"

namespace FemLib
{
             
typedef TemplateIsoparametric<FiniteElementType::LINE2, 2, FemShapeLine2, FemIntegrationGaussLine, FeExtrapolationGaussLinear> LINE2;
typedef TemplateIsoparametric<FiniteElementType::LINE3, 3, FemShapeLine3, FemIntegrationGaussLine, FeExtrapolationGaussLinear> LINE3;
typedef TemplateIsoparametric<FiniteElementType::QUAD4, 4, FemShapeQuad4, FemIntegrationGaussQuad, FeExtrapolationGaussLinear> QUAD4;
typedef TemplateIsoparametric<FiniteElementType::QUAD8, 8, FemShapeQuad8, FemIntegrationGaussQuad, FeExtrapolationGaussLinear> QUAD8;
typedef TemplateIsoparametric<FiniteElementType::QUAD9, 9, FemShapeQuad9, FemIntegrationGaussQuad, FeExtrapolationGaussLinear> QUAD9;
typedef TemplateIsoparametric<FiniteElementType::TRI3, 3, FemShapeTriangle3, FemIntegrationGaussTriangle, FeExtrapolationGaussLinear> TRI3;
typedef TemplateIsoparametric<FiniteElementType::TRI6, 6, FemShapeTriangle6, FemIntegrationGaussTriangle, FeExtrapolationGaussLinear> TRI6;

}
