
#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"


#include "TestUtil.h"

using namespace MathLib;

TEST(Math, Matrix_transposeAndMultiply)
{
    MathLib::Matrix<double> matA(1,2);
    MathLib::Matrix<double> matB(1,2);
    MathLib::Matrix<double> matC(2,2);

    matA(0,0) = 1.0;
    matA(0,1) = 2.0;
    matB(0,0) = 3.0;
    matB(0,1) = 4.0;
    matC = .0;
    matA.transposeAndMultiply(matB, matC);

    ASSERT_EQ(matC(0,0), 3.0);
    ASSERT_EQ(matC(0,1), 4.0);
    ASSERT_EQ(matC(1,0), 6.0);
    ASSERT_EQ(matC(1,1), 8.0);
}

