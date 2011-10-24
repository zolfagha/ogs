
#include <iostream>
//#include <cmath>
//#include "LinAlg/Dense/Matrix.h"
//#include "LinAlg/Dense/SymmetricMatrix.h"
//#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Sparse/SparseTableCRS.h"
//#include "LinAlg/Solvers/CG.h"
//#include "LinAlg/Dense/TemplateMatrixNd.h"
#include "LinAlg/Sparse/EigenInterface.h"
#include "Mesh.h"
#include "IO/MeshIOOGSAscii.h"
#include "IO/MeshIOLegacyVtk.h"
#include "Tools/MeshGenerator.h"
#include "CPUTimeTimer.h"
#include "RunTimeTimer.h"
#include <Eigen>
#include "MeshSparseTable.h"

#define LIS
#ifdef LIS
#include <omp.h>
#include "lis.h"
#include "LinAlg/Solvers/LisInterface.h"
#endif

using namespace std;

typedef struct {
  size_t id;
  double val;
} IndexValue;

void setDirichletBC_Case1(MeshLib::UnstructuredMesh *msh, vector<IndexValue> &list_dirichlet_bc) 
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

int main(int argc, char *argv[])
{
#ifdef LIS
    lis_initialize(&argc, &argv);
#endif
    int nthreads = 1;
    if (argc > 2)
      sscanf(argv[2], "%d", &nthreads);
    omp_set_num_threads (nthreads);
    cout << "->Start OpenMP parallelization with " << omp_get_max_threads() << " threads" << endl;
    //-- setup a problem -----------------------------------------------
    //set mesh
    std::string strMeshFile = "";
    if (argc<2) {
      std::cout << "Input a mesh file name:" << endl;
      cin >> strMeshFile;
    } else {
      strMeshFile = argv[1];
    }
    vector<MeshLib::IMesh*> vec_mesh; 
    MeshLib::MeshIOOGS::readMesh(strMeshFile, vec_mesh);
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
    MeshLib::UnstructuredMesh *msh = static_cast<MeshLib::UnstructuredMesh*>(vec_mesh.at(0));
    if (msh->getElemenet(0)->getElementType()!=MeshLib::ElementType::TRIANGLE) {
        std::cout << "Only triangle elements are supported! " << endl;
        return 0;
    }

    //set material properties
    const double K = 1.e-8;
    const double S = .0;

    //set Dirichlet BC
    vector<IndexValue> list_dirichlet_bc;
    setDirichletBC_Case1(msh, list_dirichlet_bc);

    //-- construct EQS -----------------------------------------------
    //prepare EQS
    const size_t dim_eqs = msh->getNumberOfNodes();
#define USE_EIGEN
#ifdef USE_EIGEN
    Eigen::DynamicSparseMatrix<double, Eigen::RowMajor> eqsA(dim_eqs, dim_eqs);
#else
    MathLib::SparseTableCRS<unsigned>* crs = generateSparseTableCRS<unsigned>(msh);
    //MathLib::TemplateCRSMatrix<double, int> eqsA(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
#endif

    double* eqsX(new double[dim_eqs]);
    double* eqsRHS(new double[dim_eqs]);
    for (size_t i=0; i<dim_eqs; i++) eqsX[i] = .0;
    for (size_t i=0; i<dim_eqs; i++) eqsRHS[i] = .0;


    //assembly EQS
    const size_t n_ele = msh->getNumberOfElements();
    Eigen::Matrix3d local_K = Eigen::Matrix3d::Zero();
    const size_t n_ele_nodes = 3;
    double nodes_x[n_ele_nodes], nodes_y[n_ele_nodes], nodes_z[n_ele_nodes];
    double a[n_ele_nodes], b[n_ele_nodes], c[n_ele_nodes];
    size_t dof_map[n_ele_nodes];
    MeshLib::Triangle *ele;
    double pt[3];

    cout << "->assembly EQS" << endl;

    RunTimeTimer run_timer;
    CPUTimeTimer cpu_timer;
    run_timer.start();
    cpu_timer.start();

    for (size_t i_ele=0; i_ele<n_ele; i_ele++) {
        ele = static_cast<MeshLib::Triangle*>(msh->getElemenet(i_ele));
        // setup element information
        // dof
        for (size_t i=0; i<ele->getNumberOfNodes(); i++)
            dof_map[i] = ele->getNodeID(i);
        // xyz
        for (int i=0; i<n_ele_nodes; i++) {
            msh->getNodeCoordinates(ele->getNodeID(i), pt);
            nodes_x[i] = pt[0];
            nodes_y[i] = pt[1];
            nodes_z[i] = pt[2];
        }
        // area
        const double A = 0.5*(nodes_x[0]*(nodes_y[1]-nodes_y[2])+nodes_x[1]*(nodes_y[2]-nodes_y[0])+nodes_x[2]*(nodes_y[0]-nodes_y[1]));
        // set a,b,c
        a[0] = 0.5/A*(nodes_x[1]*nodes_y[2]-nodes_x[2]*nodes_y[1]);
        b[0] = 0.5/A*(nodes_y[1]-nodes_y[2]);
        c[0] = 0.5/A*(nodes_x[2]-nodes_x[1]);
        a[1] = 0.5/A*(nodes_x[2]*nodes_y[0]-nodes_x[0]*nodes_y[2]);
        b[1] = 0.5/A*(nodes_y[2]-nodes_y[0]);
        c[1] = 0.5/A*(nodes_x[0]-nodes_x[2]);
        a[2] = 0.5/A*(nodes_x[0]*nodes_y[1]-nodes_x[1]*nodes_y[0]);
        b[2] = 0.5/A*(nodes_y[0]-nodes_y[1]);
        c[2] = 0.5/A*(nodes_x[1]-nodes_x[0]);

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
        for (int i=0; i<n_ele_nodes; i++)
            for (int j=0; j<i; j++)
                local_K(i,j) = local_K(j,i);

        // add into global EQS
        for (int i=0; i<n_ele_nodes; i++) {
            for (int j=0; j<n_ele_nodes; j++) {
#ifdef USE_EIGEN
                eqsA.coeffRef(dof_map[i], dof_map[j]) += local_K(i,j);
#else
              eqsA(dof_map[i], dof_map[j]) += local_K(i,j);
#endif
            }
        }
    }

    //outputEQS("eqs1.txt", eqsA, eqsX, eqsRHS);

    cout << "->apply BC" << endl;

    //apply Dirichlet BC
    for (size_t i=0; i<list_dirichlet_bc.size(); i++) {
        IndexValue &bc = list_dirichlet_bc.at(i);
#ifdef USE_EIGEN
        MathLib::EigenTools::setKnownXi(eqsA, eqsRHS, bc.id, bc.val);
#else
        const size_t id = bc.id;
        const double val = bc.val;
        //A(k, j) = 0.
        for (size_t j=0; j<eqsA.getNCols(); j++)
          if (eqsA(id, j)!=.0)
            eqsA(id, j) = .0;
        //A(k, k) = val,
        eqsA(id, id) = val;
        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<eqsA.getNCols(); j++)
          eqsRHS[j] -= eqsA(j, id)*val;
        //b_k = A_kk*val
        eqsRHS[id] = eqsA(id, id)*val;
        //A(i, k) = 0., i!=k
        for (size_t j=0; j<eqsA.getNCols(); j++)
          if (eqsA(j, id)!=.0 && j!=id)
            eqsA(j, id) = .0;
#endif
    }

    //apply ST
    //outputEQS("eqs2.txt", eqsA, eqsX, eqsRHS);

    //-- solve EQS -----------------------------------------------
    //setup
#ifdef USE_EIGEN
    cout << "->export Matrix for LIS" << endl;
    MathLib::SparseTableCRS<int> *crsA = MathLib::EigenTools::buildCRSMatrixFromEigenMatrix(eqsA);
#else
    MathLib::SparseTableCRS<unsigned> *crsA = crs;
#endif
    MathLib::LIS_option option;
    option.ls_method = 1;
    option.ls_precond = 0;
    option.ls_extra_arg = "";
    option.ls_max_iterations = 5000;
    option.ls_error_tolerance = 1e-10;

    //solve EQS
    cout << "->solve BC" << endl;
    CPUTimeTimer cpu_timer2;
    RunTimeTimer run_timer2;
    run_timer2.start();
    cpu_timer2.start();
    MathLib::solveWithLis(crsA, eqsX, eqsRHS, option);
    run_timer2.stop();
    cpu_timer2.stop();

    run_timer.stop();
    cpu_timer.stop();
    cout.setf(std::ios::scientific,std::ios::floatfield);
    cout.precision(12);
    cout << "== Simulation time ==" << endl;
    cout << "Total simulation:" << endl;
    cout << "CPU time = " << run_timer.elapsed() << endl;
    cout << "Run time = " << cpu_timer.elapsed() << endl;
    cout << "Linear solver:" << endl;
    cout << "CPU time = " << run_timer2.elapsed() << endl;
    cout << "Run time = " << cpu_timer2.elapsed() << endl;

    //-- output results -----------------------------------------------
    //output results
    cout << "->output results" << endl;
    std::vector<MeshLib::NodalScalarValue> nodalValues;
    string str = "Head";
    MeshLib::NodalScalarValue temp("Head", eqsX);
    nodalValues.push_back(temp);
    MeshLib::MeshIOLegacyVtk4Simulation::WriteAsciiFile("output.vtk", *msh, 1, 1.0, nodalValues);

    //release memory
#ifdef LIS
    lis_finalize();
#endif
    destroyStdVectorWithPointers(vec_mesh);
    delete crsA;
    delete [] eqsX;
    delete [] eqsRHS;

    return 0;
}

