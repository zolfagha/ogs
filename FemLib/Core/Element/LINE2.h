
#pragma once

#include "FemLib/Core/IFemElement.h"
#include "FemLib/Core/Integration.h"

namespace FemLib
{

typedef FeBaseIsoparametric<FiniteElementType::LINE2, 2, FemShapeLine2, FemIntegrationGaussLine> LINE2;

}
