
#pragma once

#include <Eigen>
#include "MathLib/LinAlg/LinearEquations/EigenLinearEquation.h"

namespace DiscreteLib
{
typedef MathLib::EigenDenseLinearEquation LocalEquation;
typedef Eigen::MatrixXd LocalMatrix;
typedef Eigen::VectorXd LocalVector;

}
