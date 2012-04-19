
#include <gtest/gtest.h>

#include <vector>

#include "MeshLib/Core/IMesh.h"
#include "GeoLib/Shape/Rectangle.h"
#include "DiscreteLib/Core/DiscreteSystem.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "FemLib/Post/Extrapolation.h"

using namespace GeoLib;
using namespace MeshLib;
using namespace DiscreteLib;
using namespace FemLib;

#include "TestUtil.h"
#include "TestExamples.h"


TEST(FEM, nonlinear1)
{

}
