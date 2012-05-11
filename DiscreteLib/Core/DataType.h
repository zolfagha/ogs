
#pragma once

#include <Eigen>
//#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/EigenLinearEquation.h"

namespace DiscreteLib
{

//typedef MathLib::Matrix<double> LocalMatrix;
typedef Eigen::MatrixXd LocalMatrix;
typedef Eigen::VectorXd LocalVector;
typedef MathLib::EigenDenseLinearEquation LocalEquation;

}
