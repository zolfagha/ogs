
#include <iostream>
//#include <cmath>
//#include "Point.h"
//#include "LinAlg/Dense/Matrix.h"
//#include "LinAlg/Dense/SymmetricMatrix.h"
//#include "LinAlg/Sparse/CRSMatrix.h"
//#include "LinAlg/Solvers/CG.h"
#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Solvers/LisInterface.h"
#include "LinAlg/Sparse/EigenInterface.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "MeshGenerator.h"
#include "CPUTimeTimer.h"
#include <Eigen>
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Sparse>

using namespace std;
using namespace MeshLib;
//using namespace MathLib;

#define TEST2

#ifdef TEST1
int main(int argc, char *argv[])
{
    cout << "Hello World!" << endl;
    StaticUnstructuredMesh<3,10,10> mesh;
    cout << "Nr. Nodes=" << mesh.getNumberOfNodes() << endl;
    cout << "Nr. Elements=" << mesh.getNumberOfElements() << endl;
    
	return 0;
}

#elif defined TEST2

int main(int argc, char *argv[])
{
    //set mesh
    std::string strMeshFile = "";
    std::cout << "Input a mesh file name:" << endl;
    cin >> strMeshFile;
    vector<IMesh*> vec_mesh; 
    MeshIOOGS::readMesh(strMeshFile, vec_mesh);
    //MeshGenerator::generateRegularMesh(2, 100, 100, 0.0, 0.0, 0.0, msh);
    //MeshIOLegacyVtk::WriteAsciiFile("test.vtk", msh);
    if (vec_mesh.size()==0) {
        std::cout << "Fail to read a mesh file: " << strMeshFile << std::endl;
        return 0;
    }

    std::cout << "== Mesh ==" << endl;
    for (size_t i=0; i<vec_mesh.size(); i++) {
        std::cout << "#Mesh " << i+1 << endl;
        std::cout << "Number of nodes   : " << vec_mesh[i]->getNumberOfNodes() << endl;
        std::cout << "Number of elements: " << vec_mesh[i]->getNumberOfElements() << endl;
    }
    std::cout << endl;
    UnstructuredMesh *msh = static_cast<UnstructuredMesh*>(vec_mesh.at(0));
    if (msh->getElemenet(0)->getElementType()!=ElementType::TRIANGLE) {
        std::cout << "Only triangle elements are supported! " << endl;
        return 0;
    }

    //set material properties
    const double K = 1.e-8;
    const double S = .0;

    //set Dirichlet BC
    typedef struct {
        size_t id;
        double val;
   } IndexValue;

    vector<IndexValue> list_dirichlet_bc;
    {
        const double head_left = 1.;
        const double head_right = 0.;
        double x_min=1e+99, x_max = -1e+99;
        double pt[3];
        //search x min/max
        for (size_t i=0; i<msh->getNumberOfNodes(); i++) {
            msh->getNodeCoordinates(i, pt);
            if (pt[0]<x_min) x_min = pt[0]; 
            if (pt[0]>x_max) x_max = pt[0]; 
        }
        //search nodes on min/max
        for (size_t i=0; i<msh->getNumberOfNodes(); i++) {
            msh->getNodeCoordinates(i, pt);
            if (abs(pt[0]-x_min)<numeric_limits<double>::epsilon()) {
                IndexValue idv={i, head_left};
                list_dirichlet_bc.push_back(idv); 
            } else if (abs(pt[0]-x_max)<numeric_limits<double>::epsilon()) {
                IndexValue idv={i, head_right};
                list_dirichlet_bc.push_back(idv); 
            }
        }
    }

    //prepare EQS
    const size_t dim_eqs = msh->getNumberOfNodes();
    Eigen::DynamicSparseMatrix<double, Eigen::RowMajor> eqsA(dim_eqs, dim_eqs);

    double* eqsX(new double[dim_eqs]);
    double* eqsRHS(new double[dim_eqs]);
    for (size_t i=0; i<dim_eqs; i++) eqsX[i] = .0;
    for (size_t i=0; i<dim_eqs; i++) eqsRHS[i] = .0;

    //-- begin solver -----------------------------------------------
    CPUTimeTimer timer;
    timer.start();

    //assembly EQS
    const size_t n_ele = msh->getNumberOfElements();
    Eigen::Matrix3d local_K;

    Triangle *ele;
    double x[3], y[3], z[3];
    double pt[3];
    double a[3], b[3], c[3];
    size_t dof_map[3];
    for (size_t i_ele=0; i_ele<n_ele; i_ele++) {
        ele = static_cast<Triangle*>(msh->getElemenet(i_ele));
        // setup element information
        // dof
        for (size_t i=0; i<ele->getNumberOfNodes(); i++)
            dof_map[i] = ele->getNodeID(i);
        // xyz
        for (int i=0; i<3; i++) {
            msh->getNodeCoordinates(ele->getNodeID(i), pt);
            x[i] = pt[0];
            y[i] = pt[1];
            z[i] = pt[2];
        }
        // area
        const double A = 0.5*(x[0]*(y[1]-y[2])+x[1]*(y[2]-y[0])+x[2]*(y[0]-y[1]));
        // set a,b,c
        a[0] = 0.5/A*(x[1]*y[2]-x[2]*y[1]);
        b[0] = 0.5/A*(y[1]-y[2]);
        c[0] = 0.5/A*(x[2]-x[1]);
        a[1] = 0.5/A*(x[2]*y[0]-x[0]*y[2]);
        b[1] = 0.5/A*(y[2]-y[0]);
        c[1] = 0.5/A*(x[0]-x[2]);
        a[2] = 0.5/A*(x[0]*y[1]-x[1]*y[0]);
        b[2] = 0.5/A*(y[0]-y[1]);
        c[2] = 0.5/A*(x[1]-x[0]);

        // assemble local EQS
        // Int{w S ph/pt + div(w) K div(p)}dA = Int{w K div(p)}dL   
        local_K(0,0) = b[0]*b[0] + c[0]*c[0];
        local_K(0,1) = b[0]*b[1] + c[0]*c[1];
        local_K(0,2) = b[0]*b[2] + c[0]*c[2];
        local_K(1,1) = b[1]*b[1] + c[1]*c[1];
        local_K(1,2) = b[1]*b[2] + c[1]*c[2];
        local_K(2,2) = b[2]*b[2] + c[2]*c[2];
        local_K *= A;
        // symmetric
        for (int i=0; i<3; i++)
            for (int j=0; j<i; j++)
                local_K(i,j) = local_K(j,i);

        // add into global EQS
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                eqsA.coeffRef(dof_map[i], dof_map[j]) += local_K(i,j);
            }
        }
    }

    //outputEQS("eqs1.txt", eqsA, eqsX, eqsRHS);

    //apply Dirichlet BC
    for (size_t i=0; i<list_dirichlet_bc.size(); i++) {
        IndexValue &bc = list_dirichlet_bc.at(i);
        MathLib::EigenTools::setKnownXi(eqsA, eqsRHS, bc.id, bc.val);
    }

    //apply ST
    //outputEQS("eqs2.txt", eqsA, eqsX, eqsRHS);

    //set CRS
    MathLib::CRS *crsA = MathLib::EigenTools::buildCRSMatrixFromEigenMatrix(eqsA);

    //solve EQS
    MathLib::LIS_option option;
    option.ls_method = 1;
    option.ls_precond = 0;
    option.ls_extra_arg = "";
    option.ls_max_iterations = 1000;
    option.ls_error_tolerance = 1e-12;
    CPUTimeTimer timer2;
    timer2.start();
    MathLib::solveWithLis(crsA, eqsX, eqsRHS, option);
    timer2.stop();

    //
    timer.stop();
    cout << "== Simulation time ==" << endl;
    cout << "Total CPU time         = " << timer.elapsed() << endl << endl;
    cout << "Linear solver CPU time = " << timer2.elapsed() << endl << endl;
    //-- end solver -------------------------------------------------

    //output results
    std::vector<MeshLib::NodalScalarValue> nodalValues;
    string str = "Head";
    MeshLib::NodalScalarValue temp("Head", eqsX);
    nodalValues.push_back(temp);
    MeshIOLegacyVtk4Simulation::WriteAsciiFile("output.vtk", *msh, 1, 1.0, nodalValues);

    //release memory
    destroyStdVectorWithPointers(vec_mesh);
    delete crsA;
    delete [] eqsX;
    delete [] eqsRHS;

    return 0;
}

#endif
