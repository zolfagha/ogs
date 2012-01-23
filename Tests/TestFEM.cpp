#include <gtest/gtest.h>

#include "FemLib/FemElement.h"
#include "FemLib/FemFunction.h"
#include "FemLib/FemFunctionProjection.h"
#include "FemLib/BC/FemDirichletBC.h"
#include "FemLib/BC/FemNeumannBC.h"
#include "FemLib/Interpolation.h"
#include "FemLib/Mapping.h"

#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Tools/MeshGenerator.h"

#include "MathLib/Vector.h"
#include "MathLib/LinearInterpolation.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/Function/Function.h"

#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/GeoGenerator.h"

#include <vector>
#include <memory>

using namespace FemLib;
using namespace GeoLib;
using namespace MeshLib;

TEST(FEM, testAll) 
{
    //#Define a problem
    //geometry
    Rectangle rec(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
    Polyline* poly_left = rec.getLeft();
    Polyline* poly_right = rec.getRight();
    //mesh
    std::auto_ptr<UnstructuredMesh2d> msh = MeshGenerator::generateRegularMesh(2, 2.0, 2, .0, .0, .0);
    //mat
    const double K = 1.e-11; // TODO: should be a function
    //discretization
    FemNodalFunctionScalar2d head(msh.get(), PolynomialOrder::Linear);
    FEMIntegrationPointFunctionVector2d vel(msh.get());
    //bc
    FemDirichletBC<double> bc1(&head, poly_right, &MathLib::FunctionConstant<double, GeoLib::Point>(.0)); //TODO should BC objects be created by fe functions?
    FemNeumannBC<double> bc2(&head, poly_left, &MathLib::FunctionConstant<double, GeoLib::Point>(1.e-5));

    //#Solve
    bc1.setup();
    bc2.setup();
    // global EQS
    const size_t n_dof = msh->getNumberOfNodes();
    MathLib::Matrix<double> globalA(n_dof, n_dof);
    std::vector<double> globalRHS(n_dof, .0);

    //assembly
    IFiniteElement *fe = head.getFiniteElement();
    MathLib::Matrix<double> localK;
    std::vector<size_t> e_node_id_list;
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement *e = msh->getElemenet(i_e);
        fe->configure(e);
        const size_t &n_dof = fe->getNumberOfDOFs();
        localK.resize(n_dof, n_dof);
        localK = .0;
        fe->computeIntTestShape(0, &localK);
        e->getNodeIDList(e_node_id_list);
        globalA.add(e_node_id_list, localK); //TODO A(id_list) += K;
    }

    //apply BC
    bc2.apply(&globalRHS);
    bc1.apply(&globalA, &globalRHS);

    //solve
    MathLib::GaussAlgorithm solver(globalA);
    std::vector<double> x(globalRHS);
    solver.execute(&x[0]); // input x contains rhs but output x contains solution.

    //update head
    head.setNodalValues(&x[0]);

    //calculate vel (vel=f(h))
    for (size_t i_e=0; i_e<msh->getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh->getElemenet(i_e);
        fe->configure(e);
        std::vector<double> local_h;
        // for each integration points
        IFemIntegration *integral = fe->getIntegrationMethod();
        IFemMapping *mapping = fe->getMapping();
        for (size_t ip=0; ip<integral->getNumberOfSamplingPoints(); ip++) {
            MathLib::Vector2D q;
            mapping->computeMappingFunctions(integral->getSamplingPoint(ip), FemMapComputation::DSHAPE);
            const MathLib::Matrix<double> *dN = mapping->getGradShapeFunction();
            dN->axpy(-K, &local_h[0], .0, q.getRawRef()); //TODO  q = - K * dN * local_h;
            vel.setIntegrationPointValue(i_e, ip, q);
        }
    }

    //for output
    FemNodalFunctionVector2d nod_vel(msh.get(), PolynomialOrder::Linear);
    mapFunctions(vel, nod_vel);
};

#if 0
void calcMass(double *pt, MathLib::Matrix<double> *mat) {
    //W^T N
    FemLagrangeElement fem; //gauss, iso, Bubnov
    double *shape = fem.computeShapeFunction(pt);
    double *test = fem.computeTestFunction(pt);
    int dof = 1;
    for (int j=0; j<dof; j++)
        for (int k=0; k<dof; k++)
            (*mat)(j,k) = test[j]*shape[k];
}

void calcLap(double *pt, MathLib::Matrix<double> *mat) {
    //dW^T dN
    FemLagrangeElement fem; //gauss, iso, Bubnov
    const MathLib::Matrix<double> *dshape = fem.computeGradShapeFunction(pt);
    MathLib::Matrix<double> *dtest = fem.computeGradTestFunction(pt);
    int dof = 1;
    //mat = (*dtest->transpose()) * (*dshape);
}

// Method 2
// - user only access FEM class
// - FEM class hides mapping, integration, test function
TEST(FEM, integral1)
{
    //input
    MeshLib::UnstructuredMesh msh;
    MeshLib::MeshGenerator::generateRegularMesh(2, 1.0, 1, .0, .0, .0, msh);
    MeshLib::IElement *e = msh.getElemenet(0);

    FemLagrangeElement fem; //gauss, iso, Bubnov
    fem.configure(e);
    size_t dof = e->getNumberOfNodes();
    //output
    MathLib::Matrix<double> M(dof, dof);
    MathLib::Matrix<double> K(dof, dof);

    //-----------------------------
    fem.integrate(calcMass, &M);
    fem.integrate(calcLap, &K);

    // M = fem.integrateDomain(W^T * S * N)
    // K = fem.integrateDomain(dW^T * K * N)
    // {Q} = fem.integrateDomain(W^T * Q)
    // {q} = fem.integrateBoundary(W^T * q)
    //
    // fdm.setEquation(M, K, Q+q); // M du/dt + K u = F
    // fdm.discretize(dt, u_i, localA, localRHS); // A = (1/dt M+ theta K), b = (1/dt M+ (1-theta) K) u_i + F 
    //
    // globalA(dof_map) += localA
    // localRHS(dof_map) += localRHS

}
#endif

