
#pragma once

#include <Eigen>
#include "MathLib/LinAlg/LinearEquations/EigenLinearEquation.h"

namespace DiscreteLib
{
typedef MathLib::EigenDenseLinearEquation LocalEquation;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;
typedef Eigen::VectorXd LocalVector;

}
