
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

#include "TestUtil.h"

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

void outputLinearEQS(MathLib::Matrix<double> &globalA, double* globalRHS)
{
    std::cout << "A=" << std::endl;
    globalA.write(std::cout);
    std::cout << "x=" << std::endl;
    for (size_t i=0; i<globalA.getNRows(); i++)
        std::cout << globalRHS[i] << " ";
    std::cout << std::endl;
}

class GWFemTestSystem
{
public:
    Rectangle *rec;
    IMesh *msh;
    MathLib::IFunction<double, double*> *K;
    FemNodalFunctionScalar *head;
    FEMIntegrationPointFunctionVector2d *vel;
    std::vector<FemDirichletBC<double>*> vec_bc1;
    std::vector<FemNeumannBC<double, double>*> vec_bc2;

    void define(MeshLib::IMesh *msh)
    {
        //#Define a problem
        //#Define a problem
        //geometry
        rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = rec->getLeft();
        Polyline* poly_right = rec->getRight();
        //mesh
        this->msh = msh;
        //discretization
        head = new FemNodalFunctionScalar(msh, PolynomialOrder::Linear);
        vel = new FEMIntegrationPointFunctionVector2d(msh);
        //bc
        vec_bc1.push_back(new FemDirichletBC<double>(head, poly_right, new MathLib::FunctionConstant<double, GeoLib::Point>(.0), new DiagonalizeMethod())); //TODO should BC objects be created by fe functions?
        vec_bc2.push_back(new FemNeumannBC<double, double>(head, poly_left, new MathLib::FunctionConstant<double, GeoLib::Point>(-1e-5)));

        K = new MathLib::FunctionConstant<double, double*>(1.e-11);
    }
    
    static void calculateHead(GWFemTestSystem &gw)
    {
        const MeshLib::IMesh *msh = gw.msh;
        for (size_t i=0; i<gw.vec_bc1.size(); i++) gw.vec_bc1[i]->setup();
        for (size_t i=0; i<gw.vec_bc2.size(); i++) gw.vec_bc2[i]->setup();
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

        outputLinearEQS(*globalA, globalRHS);

        //apply BC
        for (size_t i=0; i<gw.vec_bc2.size(); i++) gw.vec_bc2[i]->apply(globalRHS);
        //outputLinearEQS(globalA, globalRHS);
        for (size_t i=0; i<gw.vec_bc1.size(); i++) gw.vec_bc1[i]->apply(eqs);
        //outputLinearEQS(*globalA, globalRHS);

        //solve
        eqs.solve();

        //update head
        gw.head->setNodalValues(eqs.getX());
    }

    static void calculateVelocity(GWFemTestSystem &gw)
    {
        const MeshLib::IMesh *msh = gw.msh;
        //calculate vel (vel=f(h))
        for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
            MeshLib::IElement* e = msh->getElemenet(i_e);
            IFiniteElement *fe = gw.head->getFiniteElement(e);
            std::vector<double> local_h(e->getNumberOfNodes());
            for (size_t j=0; j<e->getNumberOfNodes(); j++)
                local_h[j] = gw.head->getValue(e->getNodeID(j));
            // for each integration points
            IFemNumericalIntegration *integral = fe->getIntegrationMethod();
            double r[2] = {};
            const size_t n_gp = integral->getNumberOfSamplingPoints();
            gw.vel->setNumberOfIntegationPoints(i_e, n_gp);
            std::vector<double> xi(e->getNumberOfNodes());
            std::vector<double> yi(e->getNumberOfNodes());
            for (size_t i=0; i<e->getNumberOfNodes(); i++) {
                GeoLib::Point &pt = msh->getNodeCoordinates(e->getNodeID(i));
                xi[i] = pt[0];
                yi[i] = pt[1];
            }
            for (size_t ip=0; ip<n_gp; ip++) {
                MathLib::Vector2D q;
                q.getRawRef()[0] = .0;
                q.getRawRef()[1] = .0;
                integral->getSamplingPoint(ip, r);
                fe->computeBasisFunctions(r);
                const MathLib::Matrix<double> *dN = fe->getGradBasisFunction();
                MathLib::Matrix<double>*N = fe->getBasisFunction();
                std::vector<double> xx(2); 
                N->axpy(1.0, &xi[0], .0, &xx[0]);
                N->axpy(1.0, &yi[0], .0, &xx[1]);

                double k = gw.K->eval(&xx[0]);
                dN->axpy(-k, &local_h[0], .0, q.getRawRef()); //TODO  q = - K * dN * local_h;
                gw.vel->setIntegrationPointValue(i_e, ip, q);
            }
        }
    }
};





void getGWExpectedHead(std::vector<double> &expected)
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
    GWFemTestSystem gw;
    MeshLib::IMesh *msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTestSystem::calculateHead(gw);

    double *h = gw.head->getNodalValues();
    std::vector<double> expected;
    getGWExpectedHead(expected);

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h, gw.head->getNumberOfNodes());
}

TEST(FEM, testStructuredMesh)
{
    GWFemTestSystem gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    //#Solve
    GWFemTestSystem::calculateHead(gw);

    double *h = gw.head->getNodalValues();
    std::vector<double> expected;
    getGWExpectedHead(expected);

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h, gw.head->getNumberOfNodes());
}

TEST(FEM, ExtrapolateAverage1)
{
    GWFemTestSystem gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    GWFemTestSystem::calculateHead(gw);
    GWFemTestSystem::calculateVelocity(gw);

    FemNodalFunctionVector2d nodal_vel(gw.msh, PolynomialOrder::Linear);
    FemExtrapolationAverage<MathLib::Vector2D> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);

    MathLib::Vector2D *v = nodal_vel.getNodalValues();
    ASSERT_DOUBLE_ARRAY_EQ(MathLib::Vector2D(1.e-5, .0), v, gw.head->getNumberOfNodes());
}

template<typename Tval, typename Tpos>
class MyFunction : public MathLib::IFunction<Tval, Tpos>
{
public:
    virtual Tval eval(const Tpos& x) 
    {
        if (x[0]<1.0) return 1e-11;
        else return 2e-11;
    };

    MathLib::IFunction<Tval, Tpos>* clone() const {return 0;};

};

TEST(FEM, ExtrapolateAverage2)
{
    GWFemTestSystem gw;
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    gw.define(msh);
    delete gw.K;
    gw.K = new MyFunction<double, double*>();
    gw.vec_bc1.push_back(new FemDirichletBC<double>(gw.head, gw.rec->getLeft(), new MathLib::FunctionConstant<double, GeoLib::Point>(2.e+6), new DiagonalizeMethod())); 
    delete gw.vec_bc2[0];
    gw.vec_bc2.clear();

    GWFemTestSystem::calculateHead(gw);
    GWFemTestSystem::calculateVelocity(gw);

    FemNodalFunctionVector2d nodal_vel(gw.msh, PolynomialOrder::Linear);
    FemExtrapolationAverage<MathLib::Vector2D> extrapo;
    extrapo.extrapolate(*gw.vel, nodal_vel);

    std::vector<double> exH(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) exH[i] = 2.e+6;
        if (i%3==1) exH[i] = 2./3.*1.e+6;
        if (i%3==2) exH[i] = 0.e+6;
    }

    double *h = gw.head->getNodalValues();
    ASSERT_DOUBLE_ARRAY_EQ(&exH[0], h, gw.head->getNumberOfNodes());

    MathLib::Vector2D *v = nodal_vel.getNodalValues();
    ASSERT_DOUBLE_ARRAY_EQ(MathLib::Vector2D(4./3.*1.e-5, .0), v, gw.head->getNumberOfNodes());
}

