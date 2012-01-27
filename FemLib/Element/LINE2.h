
#pragma once

#include "FemLib/IFemElement.h"
#include "FemLib/Integration.h"

namespace FemLib
{

typedef FeBaseIsoparametric<FiniteElementType::LINE2, 2, FemShapeLine2, FemIntegrationGaussLine> LINE2;

}
