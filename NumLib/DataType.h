
#pragma once

#include <Eigen>
#include "MathLib/LinAlg/LinearEquations/EigenLinearEquation.h"

namespace NumLib
{
typedef Eigen::MatrixXd LocalMatrix;
typedef Eigen::VectorXd LocalVector;
typedef MathLib::EigenDenseLinearEquation LocalEquation;
} //end
