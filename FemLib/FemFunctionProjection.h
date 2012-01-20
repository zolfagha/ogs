
#pragma once

#include "FemFunction.h"
#include "Extrapolation.h"

namespace FemLib
{

template<typename Tvalue, typename Tpos>
void mapFunctions(FEMIntegrationPointFunction<Tvalue, Tpos> &ele, FEMNodalFunction<Tvalue, Tpos> &nod)
{
    IFemExtrapolation<Tvalue, Tpos> *method = FEMExtrapolationFactory<Tvalue, Tpos>::create(FEMExtrapolationMethod::Linear);
    method->extrapolate(ele, nod);
};


}
