
#include "Fem.h"
#include "FemFunction.h"
#include "BoundaryConditions.h"

#include "MeshLib/Core/Mesh.h"
#include "MeshLib/Tools/MeshGenerator.h"

#include "MathLib/Vector.h"
#include "MathLib/LinearInterpolation.h"

#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

#include <vector>

using namespace FemLib;

void testFunction() 
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

    //define discretization
    FEMNodalFunction<double, MathLib::Vector2D> head(&msh, LagrangeOrder::Linear);
    FEMIntegrationPointFunction<MathLib::Vector2D, MathLib::Vector2D> vel(&msh);

    //prepare bc
    DirichletBC bc;
    bc.set(&head, &poly_left, &MathLib::DistributionConstant(.0));
    NeumannBC st;

    // global EQS
    const size_t n_dof = msh.getNumberOfNodes();
    MathLib::Matrix<double> globalA(n_dof, n_dof);
    std::vector<double> globalRHS(n_dof, .0);

    //assembly
    IFiniteElement *fe = head.getFiniteElement();
    MathLib::Matrix<double> localK;
    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        const size_t e_nnodes = e->getNumberOfNodes();
        localK.resize(e_nnodes, e_nnodes);
        localK = .0;
        fe->configure(e);
        fe->integrateShapeShape(0, &localK);
    }

    //apply BC

    //solve

    //calculate vel
    //map gauss points to nodal values

}

void testIntegration1() 
{

    //input
    MeshLib::IElement *e;
    FemLagrangeElement fem; //gauss, iso, Bubnov
    fem.configure(e);
    size_t dof = 3;
    //output
    MathLib::Matrix<double> M(dof, dof);
    MathLib::Matrix<double> K(dof, dof);

    //-----------------------------
    // Int{N(x)^T N(x)}dV
    // = Int{N(s)^T N(s) j(s)}dS
    // = sum_i {N^T(si)N(si)j(si)w(si)}
    FemIntegrationGauss *gauss = static_cast<FemIntegrationGauss*>(fem.getIntegration());
    for (int i=0; i<gauss->getNumberOfSamplingPoints(); i++) {
        //fem.getShapeFunctions(i, shape, dshape, &j_det);
        double *shape = fem.getShapeFunction(i);
        MathLib::Matrix<double> *dshape = fem.getGradShapeFunction(i);
        double *test = fem.getTestFunction(i);
        MathLib::Matrix<double> *dtest = fem.getGradTestFunction(i);
        double fkt = gauss->getWeight(i)*fem.getDetJacobian(i);

        //M
        for (int j=0; j<dof; j++)
            for (int k=0; k<dof; k++)
                M(j,k) += test[j]*shape[k]*fkt;

        //K
        //K+= dshape^T dtest * fkt
    }
    
}

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
    MathLib::Matrix<double> *dshape = fem.computeGradShapeFunction(pt);
    MathLib::Matrix<double> *dtest = fem.computeGradTestFunction(pt);
    int dof = 1;
    //mat = (*dtest->transpose()) * (*dshape);
}

// Method 2
// - user only access FEM class
// - FEM class hides mapping, integration, test function
void testIntegration2()
{
    //input
    MeshLib::IElement *e;
    FemLagrangeElement fem; //gauss, iso, Bubnov
    fem.configure(e);
    size_t dof = 3;
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


int main(int argc, char* argv[]) {

    testIntegration1();
    testIntegration2();
    testFunction();

    return 0;
}
