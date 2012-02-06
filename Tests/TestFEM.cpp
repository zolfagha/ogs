
#include <gtest/gtest.h>

#include "MathLib/Vector.h"
#include "MathLib/LinearInterpolation.h"
//#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/Function/FemFunctionProjection.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "FemLib/Post/Extrapolation.h"

#include <vector>
#include <memory>

using namespace FemLib;
using namespace GeoLib;
using namespace MeshLib;

typedef MathLib::Matrix<double> GlobalMatrixType;
typedef std::vector<double> GlobalVectorType;

void outputLinearEQS(MathLib::Matrix<double> &globalA, std::vector<double> &globalRHS)
{
    std::cout << "A=" << std::endl;
    globalA.write(std::cout);
    std::cout << "x=" << std::endl;
    for (size_t i=0; i<globalRHS.size(); i++)
        std::cout << globalRHS[i] << " ";
    std::cout << std::endl;
}

class GWProblem
{
public:
    IMesh *msh;
    MathLib::IFunction<double, double*> *K;
    FemNodalFunctionScalar *head;
    FEMIntegrationPointFunctionVector2d *vel;
    FemDirichletBC<double> *bc1;
    FemNeumannBC<double> *bc2;
};

void solveGW(GWProblem &gw)
{
    const MeshLib::IMesh *msh = gw.msh;
    gw.bc1->setup();
    gw.bc2->setup();
    // global EQS
    const size_t n_dof = msh->getNumberOfNodes();
    MathLib::DenseLinearEquations eqs;
    eqs.create(n_dof);

    MathLib::DenseLinearEquations::MatrixType* globalA = eqs.getA();
    double* globalRHS = eqs.getRHS();

    //assembly
    MathLib::Matrix<double> localK;
    std::vector<size_t> e_node_id_list;
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement *e = msh->getElemenet(i_e);
        IFiniteElement *fe = gw.head->getFiniteElement(e);
        const size_t &n_dof = fe->getNumberOfVariables();
        localK.resize(n_dof, n_dof);
        localK = .0;
        fe->integrateDWxDN(gw.K, localK);
        e->getNodeIDList(e_node_id_list);
        globalA->add(e_node_id_list, localK); //TODO A(id_list) += K;
    }

    //outputLinearEQS(globalA, globalRHS);

    //apply BC
    gw.bc2->apply(globalRHS);
    //outputLinearEQS(globalA, globalRHS);
    gw.bc1->apply(eqs);
    //outputLinearEQS(globalA, globalRHS);

    //solve
    eqs.solve();

    //update head
    gw.head->setNodalValues(eqs.getX());

    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh->getElemenet(i_e);
        IFiniteElement *fe = gw.head->getFiniteElement(e);
        std::vector<double> local_h(e->getNumberOfNodes());
        for (size_t j=0; j<e->getNumberOfNodes(); j++)
            local_h[j] = gw.head->getValue(e->getNodeID(j));
        // for each integration points
        IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        double x[2] = {};
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        gw.vel->setNumberOfIntegationPoints(i_e, n_gp);
        for (size_t ip=0; ip<n_gp; ip++) {
            MathLib::Vector2D q;
            q.getRawRef()[0] = .0;
            q.getRawRef()[1] = .0;
            integral->getSamplingPoint(ip, x);
            fe->computeBasisFunctions(x);
            const MathLib::Matrix<double> *dN = fe->getGradBasisFunction();
            double k = gw.K->eval(x);
            dN->axpy(-k, &local_h[0], .0, q.getRawRef()); //TODO  q = - K * dN * local_h;
            gw.vel->setIntegrationPointValue(i_e, ip, q);
        }
    }
}

TEST(FEM, testUnstructuredMesh)
{
    //#Define a problem
    GWProblem gw;
    //#Define a problem
    //geometry
    Rectangle rec(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
    Polyline* poly_left = rec.getLeft();
    Polyline* poly_right = rec.getRight();
    //mesh
    gw.msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    //discretization
    gw.head = new FemNodalFunctionScalar(gw.msh, PolynomialOrder::Linear);
    gw.vel = new FEMIntegrationPointFunctionVector2d(gw.msh);
    //bc
    gw.bc1 = new FemDirichletBC<double>(gw.head, poly_right, &MathLib::FunctionConstant<double, GeoLib::Point>(.0), &DiagonalizeMethod()); //TODO should BC objects be created by fe functions?
    gw.bc2 = new FemNeumannBC<double>(gw.head, poly_left, &MathLib::FunctionConstant<double, GeoLib::Point>(1.e-5));

    gw.K = new MathLib::FunctionConstant<double, double*>(1.e-11);

    //#Solve
    solveGW(gw);
}

TEST(FEM, testStructuredMesh)
{
    //#Define a problem
    GWProblem gw;
    //#Define a problem
    //geometry
    Rectangle rec(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
    Polyline* poly_left = rec.getLeft();
    Polyline* poly_right = rec.getRight();
    //mesh
    gw.msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    //discretization
    gw.head = new FemNodalFunctionScalar(gw.msh, PolynomialOrder::Linear);
    gw.vel = new FEMIntegrationPointFunctionVector2d(gw.msh);
    //bc
    gw.bc1 = new FemDirichletBC<double>(gw.head, poly_right, &MathLib::FunctionConstant<double, GeoLib::Point>(.0), &DiagonalizeMethod()); //TODO should BC objects be created by fe functions?
    gw.bc2 = new FemNeumannBC<double>(gw.head, poly_left, &MathLib::FunctionConstant<double, GeoLib::Point>(1.e-5));

    gw.K = new MathLib::FunctionConstant<double, double*>(1.e-11);

    //#Solve
    solveGW(gw);

}

TEST(FEM, ExtrapolateAverage)
{
    //#Define a problem
    GWProblem gw;
    //#Define a problem
    //geometry
    Rectangle rec(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
    Polyline* poly_left = rec.getLeft();
    Polyline* poly_right = rec.getRight();
    //mesh
    gw.msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    //discretization
    gw.head = new FemNodalFunctionScalar(gw.msh, PolynomialOrder::Linear);
    gw.vel = new FEMIntegrationPointFunctionVector2d(gw.msh);
    //bc
    gw.bc1 = new FemDirichletBC<double>(gw.head, poly_right, &MathLib::FunctionConstant<double, GeoLib::Point>(.0), &DiagonalizeMethod()); //TODO should BC objects be created by fe functions?
    gw.bc2 = new FemNeumannBC<double>(gw.head, poly_left, &MathLib::FunctionConstant<double, GeoLib::Point>(1.e-5));

    gw.K = new MathLib::FunctionConstant<double, double*>(1.e-11);
    solveGW(gw);

    FemNodalFunctionVector2d nodal_vel(gw.msh, PolynomialOrder::Linear);
    FEMExtrapolationAverage<MathLib::Vector2D> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);
}

TEST(FEM, ExtrapolateLinear)
{

}
