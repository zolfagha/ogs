
#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

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


TEST(Math, MatrixAddSubMatrix)
{
    Matrix<double> _m(4,4);
    _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
    _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
    _m(2,2) = 4.0; _m(2,3) = -1.0;
    _m(3,3) = 4.0;
    for (size_t i=0; i<4; i++)
        for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);

    DenseLinearEquations eqs;
    eqs.create(10);

    for (size_t i=0; i<3; i++) {
        std::vector<long> pos(4);
        for (size_t j=0; j<4; j++)
            pos[j] = i*3 + j;
        eqs.addAsub(pos, _m);
    }

    double expected_m[] = {
         4, -1, -2, -1,  0,  0,  0,  0,  0,  0,
        -1,  4, -1, -2,  0,  0,  0,  0,  0,  0,
        -2, -1,  4, -1,  0,  0,  0,  0,  0,  0,
        -1, -2, -1,  8, -1, -2, -1,  0,  0,  0,
         0,  0,  0, -1,  4, -1, -2,  0,  0,  0,
         0,  0,  0, -2, -1,  4, -1,  0,  0,  0,
         0,  0,  0, -1, -2, -1,  8, -1, -2, -1,
         0,  0,  0,  0,  0,  0, -1,  4, -1, -2,
         0,  0,  0,  0,  0,  0, -2, -1,  4, -1,
         0,  0,  0,  0,  0,  0, -1, -2, -1,  4
    };

    ASSERT_DOUBLE_ARRAY_EQ(expected_m, eqs.getA()->getData(), 100);
}