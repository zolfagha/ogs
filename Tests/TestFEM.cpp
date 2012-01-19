#include <gtest/gtest.h>

#include "FemLib/Fem.h"
#include "FemLib/FemFunction.h"
#include "FemLib/BoundaryConditions.h"
#include "FemLib/Projection.h"

#include "MeshLib/Core/Mesh.h"
#include "MeshLib/Tools/MeshGenerator.h"

#include "MathLib/Vector.h"
#include "MathLib/LinearInterpolation.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"

#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

#include <vector>

using namespace FemLib;

int add (int x, int y) {return x+y;};

TEST(AddTest, Test1)
{
    ASSERT_EQ(2, add(1, 1));
}

TEST(FEM, testAll) 
{
    //geo
    std::vector<GeoLib::Point*> pnt_vec;
    pnt_vec.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
    pnt_vec.push_back(new GeoLib::Point(2.0, 0.0, 0.0));
    pnt_vec.push_back(new GeoLib::Point(2.0, 2.0, 0.0));
    pnt_vec.push_back(new GeoLib::Point(0.0, 2.0, 0.0));
    GeoLib::Polyline poly_left(pnt_vec);
    poly_left.addPoint(0);
    poly_left.addPoint(3);
    GeoLib::Polyline poly_right(pnt_vec);
    poly_left.addPoint(1);
    poly_left.addPoint(2);
    //mesh
    MeshLib::UnstructuredMesh msh;
    MeshLib::MeshGenerator::generateRegularMesh(2, 2.0, 2, .0, .0, .0, msh);
    //mat
    const double K = 1.e-11;

    //define discretization
    FEMNodalFunction<double, MathLib::Vector2D> head(&msh, LagrangeOrder::Linear);
    FEMIntegrationPointFunction<MathLib::Vector2D, MathLib::Vector2D> vel(&msh);

    //prepare bc
    DirichletBC bc1;
    bc1.set(&head, &poly_right, &MathLib::DistributionConstant(.0));
    NeumannBC bc2;
    bc2.set(&head, &poly_left, &MathLib::DistributionConstant(1.e-5));

    // global EQS
    const size_t n_dof = msh.getNumberOfNodes();
    MathLib::Matrix<double> globalA(n_dof, n_dof);
    std::vector<double> globalRHS(n_dof, .0);

    //assembly
    IFiniteElement *fe = head.getFiniteElement();
    MathLib::Matrix<double> localK;
    std::vector<size_t> e_node_id_list;
    for (size_t i_e=0; i_e<msh.getNumberOfElements(); i_e++) {
        MeshLib::IElement *e = msh.getElemenet(i_e);
        const size_t e_nnodes = e->getNumberOfNodes();
        localK.resize(e_nnodes, e_nnodes);
        localK = .0;
        fe->configure(e);
        fe->integrateShapeShape(0, &localK);
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
    for (size_t i_e=0; i_e<msh.getNumberOfElements(); i_e++) {
        MeshLib::IElement* e = msh.getElemenet(i_e);
        fe->configure(e);
        std::vector<double> local_h;
        // for each integration points
        IFemIntegration *integral = fe->getIntegrationMethod();
        for (size_t ip=0; ip<integral->getNumberOfSamplingPoints(); ip++) {
            MathLib::Vector2D q;
            const MathLib::Matrix<double> *dN = fe->computeGradShapeFunction(integral->getSamplingPoint(ip));
            dN->axpy(-K, &local_h[0], .0, q.getRaw()); //TODO  q = - K * dN * local_h;
            vel.setIntegrationPointValue(i_e, ip, q);
        }
    }

    //for output
    FEMNodalFunction<MathLib::Vector2D, MathLib::Vector2D> nod_vel(&msh, LagrangeOrder::Linear);
    mapFunctions(vel, nod_vel);
};

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

