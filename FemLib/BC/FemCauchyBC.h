
#pragma once

#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemFunction.h"
#include "IFemBC.h"

namespace FemLib
{
class CauchyBC : IFemBC
{
public:
    void set(TemplateFEMNodalFunction<double,double> *fem, int geo, int func);
};

}
