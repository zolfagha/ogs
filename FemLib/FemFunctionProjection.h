
#pragma once

#include "FemFunction.h"
#include "Extrapolation.h"

namespace FemLib
{

template<typename Tvalue>
void mapFunctions(TemplateFEMIntegrationPointFunction<Tvalue> &ele, TemplateFEMNodalFunction<Tvalue> &nod)
{
    IFemExtrapolation<Tvalue> *method = FEMExtrapolationFactory<Tvalue>::create(FEMExtrapolationMethod::Linear);
    method->extrapolate(ele, nod);
};


}
