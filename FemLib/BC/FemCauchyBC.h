
#pragma once

#include "FemFunction.h"

#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

namespace FemLib
{
class CauchyBC
{
public:
    void set(FEMNodalFunction<double,double> *fem, int geo, int func);
};

}
