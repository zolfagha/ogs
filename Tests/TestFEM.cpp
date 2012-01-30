
#include <gtest/gtest.h>

#include "MathLib/Vector.h"
#include "MathLib/LinearInterpolation.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/Function/Function.h"

#include "GeoLib/Core/Point.h"
#include "GeoLib/Core/Polyline.h"
#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Tools/MeshGenerator.h"

#include "FemLib/Core/IFemElement.h"
#include "FemLib/Function/FemFunction.h"
#include "FemLib/Function/FemFunctionProjection.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "FemLib/Core/FemCoorinatesMapping.h"

#include "ModelLib/GROUNDWATER_FLOW.h"

#include <vector>
#include <memory>

using namespace FemLib;
using namespace GeoLib;
using namespace MeshLib;

void outputLinearEQS(MathLib::Matrix<double> &globalA, std::vector<double> &globalRHS)
{
    std::cout << "A=" << std::endl;
    globalA.write(std::cout);
    std::cout << "x=" << std::endl;
    for (size_t i=0; i<globalRHS.size(); i++)
        std::cout << globalRHS[i] << " ";
    std::cout << std::endl;
}

TEST(MODEL, test1)
{
    ModelLib::GROUNDWATER_FLOW gw;
}

TEST(FEM, testAll) 
{
    typedef MathLib::Matrix<double> GlobalMatrixType;
    typedef std::vector<double> GlobalVectorType;

    //#Define a problem
    //geometry
    Rectangle rec(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
    Polyline* poly_left = rec.getLeft();
    Polyline* poly_right = rec.getRight();
    //mesh
    std::auto_ptr<UnstructuredMesh2d> msh = MeshGenerator::generateRegularQuadMesh(2.0, 2, .0, .0, .0);
    //mat
    const double K = 1.e-11; // TODO: should be a function
    //discretization
    FemNodalFunctionScalar2d head(msh.get(), PolynomialOrder::Linear);
    FEMIntegrationPointFunctionVector2d vel(msh.get());
    //bc
    FemDirichletBC<double> bc1(&head, poly_right, &MathLib::FunctionConstant<double, GeoLib::Point>(.0), &DiagonalizeMethod<GlobalMatrixType,GlobalVectorType,double>()); //TODO should BC objects be created by fe functions?
    FemNeumannBC<double> bc2(&head, poly_left, &MathLib::FunctionConstant<double, GeoLib::Point>(1.e-5));

    //#Solve
    bc1.setup();
    bc2.setup();
    // global EQS
    const size_t n_dof = msh->getNumberOfNodes();
    GlobalMatrixType globalA(n_dof, n_dof, .0);
    GlobalVectorType globalRHS(n_dof, .0);

    //assembly
    MathLib::Matrix<double> localK;
    std::vector<size_t> e_node_id_list;
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement *e = msh->getElemenet(i_e);
        IFiniteElement *fe = head.getFiniteElement(e);
        const size_t &n_dof = fe->getNumberOfVariables();
        localK.resize(n_dof, n_dof);
        localK = .0;
        fe->integrateDWxDN(&MathLib::FunctionConstant<double, double*>(K), localK);
        e->getNodeIDList(e_node_id_list);
        globalA.add(e_node_id_list, localK); //TODO A(id_list) += K;
    }

    //outputLinearEQS(globalA, globalRHS);

    //apply BC
    bc2.apply(&globalRHS);
    //outputLinearEQS(globalA, globalRHS);
    bc1.apply(&globalA, &globalRHS);
    //outputLinearEQS(globalA, globalRHS);

    //solve
    MathLib::GaussAlgorithm solver(globalA);
    GlobalVectorType x(globalRHS);
    solver.execute(&x[0]); // input x contains rhs but output x contains solution.

    //update head
    head.setNodalValues(&x[0]);

    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh->getElemenet(i_e);
        IFiniteElement *fe = head.getFiniteElement(e);
        std::vector<double> local_h(e->getNumberOfNodes());
        for (size_t j=0; j<e->getNumberOfNodes(); j++)
            local_h[j] = head.getValue(e->getNodeID(j));
        // for each integration points
        IFemNumericalIntegration *integral = fe->getIntegrationMethod();
        double x[2] = {};
        const size_t n_gp = integral->getNumberOfSamplingPoints();
        vel.setNumberOfIntegationPoints(i_e, n_gp);
        for (size_t ip=0; ip<n_gp; ip++) {
            MathLib::Vector2D q;
            q.getRawRef()[0] = .0;
            q.getRawRef()[1] = .0;
            integral->getSamplingPoint(ip, x);
            fe->computeBasisFunctions(x);
            const MathLib::Matrix<double> *dN = fe->getGradBasisFunction();
            dN->axpy(-K, &local_h[0], .0, q.getRawRef()); //TODO  q = - K * dN * local_h;
            vel.setIntegrationPointValue(i_e, ip, q);
        }
    }

};

