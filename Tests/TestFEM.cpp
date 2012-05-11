
#include <gtest/gtest.h>

#include "MathLib/Vector.h"
#include "MathLib/Interpolation/LinearInterpolation.h"
//#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "NumLib/Function/Function.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"

#include "DiscreteLib/Core/DiscreteSystem.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "FemLib/Post/Extrapolation.h"

#include "TestExamples.h"
#include "TestUtil.h"

#include <vector>
#include <memory>

using namespace FemLib;
using namespace GeoLib;
using namespace MeshLib;
using namespace DiscreteLib;

typedef MathLib::Matrix<double> GlobalMatrixType;
typedef std::vector<double> GlobalVectorType;

//static void outputLinearEQS(MathLib::Matrix<double> &globalA, std::vector<double> &globalRHS)
//{
//    std::cout << "A=" << std::endl;
//    globalA.write(std::cout);
//    std::cout << "x=" << std::endl;
//    for (size_t i=0; i<globalRHS.size(); i++)
//        std::cout << globalRHS[i] << " ";
//    std::cout << std::endl;
//}
//
//static void outputLinearEQS(MathLib::Matrix<double> &globalA, double* globalRHS)
//{
//    std::cout << "A=" << std::endl;
//    globalA.write(std::cout);
//    std::cout << "x=" << std::endl;
//    for (size_t i=0; i<globalA.getNRows(); i++)
//        std::cout << globalRHS[i] << " ";
//    std::cout << std::endl;
//}




static void getGWExpectedHead(std::vector<double> &expected)
{
    expected.resize(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }
}

TEST(FEM, testUnstructuredMesh)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTest::calculateHead(gw);

    DiscreteLib::DiscreteVector<double>* h = gw.head->getNodalValues();
    std::vector<double> expected;
    getGWExpectedHead(expected);

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h)[0], gw.head->getNumberOfNodes());
}

TEST(FEM, testStructuredMesh)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTest::calculateHead(gw);

    DiscreteLib::DiscreteVector<double>* h = gw.head->getNodalValues();
    std::vector<double> expected;
    getGWExpectedHead(expected);

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], &(*h)[0], gw.head->getNumberOfNodes());
}



TEST(FEM, ExtrapolateAverage1)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    GWFemTest::calculateHead(gw);
    GWFemTest::calculateVelocity(gw);

    FemNodalFunctionVector2d nodal_vel(*gw.dis, PolynomialOrder::Linear);
    FemExtrapolationAverage<NumLib::LocalVector> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);

    DiscreteLib::DiscreteVector<NumLib::LocalVector> *v = nodal_vel.getNodalValues();
    NumLib::LocalVector expected(2);
    expected[0] = 1.e-5;
    expected[1] = .0;
    ASSERT_DOUBLE_ARRAY_EQ(expected, (*v)[0], expected.size());
}

class MyFunction : public NumLib::ITXFunction
{
public:
	virtual ~MyFunction() {};
    virtual void eval(const double* x, double &v)
    {
        if (x[0]<1.0) v = 1e-11;
        else v = 2e-11;
    };

    MyFunction* clone() const {return 0;};

};

TEST(FEM, ExtrapolateAverage2)
{
    GWFemTest gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    delete gw._K;
    gw._K = new MyFunction();
    gw.vec_bc1.push_back(new FemDirichletBC(msh, gw.rec->getLeft(), new NumLib::TXFunctionConstant(2.e+6)));
    delete gw.vec_bc2[0];
    gw.vec_bc2.clear();

    GWFemTest::calculateHead(gw);
    GWFemTest::calculateVelocity(gw);

    FemNodalFunctionVector2d nodal_vel(*gw.dis, PolynomialOrder::Linear);
    FemExtrapolationAverage<NumLib::LocalVector> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);

    std::vector<double> exH(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) exH[i] = 2.e+6;
        if (i%3==1) exH[i] = 2./3.*1.e+6;
        if (i%3==2) exH[i] = 0.e+6;
    }

    DiscreteLib::DiscreteVector<double> *h = gw.head->getNodalValues();
    ASSERT_DOUBLE_ARRAY_EQ(&exH[0], &(*h)[0], gw.head->getNumberOfNodes());

    DiscreteLib::DiscreteVector<NumLib::LocalVector> *v = nodal_vel.getNodalValues();
    NumLib::LocalVector expected(2);
    expected[0] = 4./3.*1.e-5;
    expected[1] = .0;
    ASSERT_DOUBLE_ARRAY_EQ(expected, (*v)[0], expected.size());
}

