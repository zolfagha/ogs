
#include <iostream>
//#include <cmath>
//#include "LinAlg/Dense/Matrix.h"
//#include "LinAlg/Dense/SymmetricMatrix.h"
#include "LinAlg/Sparse/CRSMatrixOpenMP.h"
#include "LinAlg/Sparse/CRSMatrixDiagPrecond.h"
#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Sparse/SparseTableCRS.h"

#include "sparse.h"

//#include "LinAlg/Dense/TemplateMatrixNd.h"
#include "LinAlg/Sparse/EigenInterface.h"
#include "LinAlg/Solvers/CG.h"
#include "Mesh.h"
#include "IO/MeshIOOGSAscii.h"
#include "IO/MeshIOLegacyVtk.h"
#include "Tools/MeshGenerator.h"
#include "CPUTimeTimer.h"
#include "RunTimeTimer.h"
#include <Eigen>
#include "MeshSparseTable.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

//-----------------------------------------------------------------------------
// Configuration
//-----------------------------------------------------------------------------
//#define LIS
//#define USE_EIGEN
#define CRS_MATRIX
//#define NO_SOLVER
//#define OUTPUT_VTK

#ifdef LIS
#include "lis.h"
#include "LinAlg/Solvers/LisInterface.h"
#define INDEX_TYPE int
#else
#define INDEX_TYPE unsigned
#endif

//-----------------------------------------------------------------------------
// Tools for testing purpose
//-----------------------------------------------------------------------------
struct IndexValue {
	IndexValue(size_t idx, double value) :
		id(idx), val(value)
	{};

	size_t id;
	double val;
};

void setDirichletBC_Case1(MeshLib::UnstructuredMesh *msh, vector<IndexValue> &list_dirichlet_bc)
{
  const double head_left = 1.0;
  const double head_right = 0.0;
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
      list_dirichlet_bc.push_back(IndexValue(i, head_left));
    } else if (abs(pt[0]-x_max)<numeric_limits<double>::epsilon()) {
      list_dirichlet_bc.push_back(IndexValue(i, head_right));
    }
  }
}

#ifndef USE_EIGEN
void setKnownXi_ReduceSizeOfEQS(vector<IndexValue> &list_dirichlet_bc, MathLib::CRSMatrix<double, INDEX_TYPE> &eqsA, double* org_eqsRHS, double* org_eqsX, double** eqsRHS, double** eqsX, map<INDEX_TYPE,INDEX_TYPE> &map_solved_orgEqs)
{
	std::cout << "first step of applying BC ... " << std::flush;
	RunTimeTimer run;
	run.start();
    const size_t n_org_rows = eqsA.getNRows();
    vector<INDEX_TYPE> removed_rows(list_dirichlet_bc.size());
    for (size_t i=0; i<list_dirichlet_bc.size(); i++) {
        IndexValue &bc = list_dirichlet_bc.at(i);
        const size_t id = bc.id;
        const double val = bc.val;
        removed_rows.at(i) = id;

        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<eqsA.getNCols(); j++)
            org_eqsRHS[j] -= eqsA.getValue(j, id)*val;
        //b_k = A_kk*val
        org_eqsRHS[id] = val; //=eqsA(id, id)*val;
        org_eqsX[id] = val; //=eqsA(id, id)*val;
    }
    run.stop();
    std::cout << run.elapsed() << " s" << std::endl;

    std::cout << "second step of applying BC (erase entries from matrix) ... " << std::flush;
	run.start();
    //remove rows and columns
    eqsA.eraseEntries(removed_rows.size(), &removed_rows[0]);
    run.stop();
    std::cout << run.elapsed() << " s" << std::endl;

    std::cout << "third step of applying BC (create new x and RHS) ... " << std::flush;
	run.start();
    //remove X,RHS
    (*eqsX) = new double[n_org_rows-removed_rows.size()];
    (*eqsRHS) = new double[n_org_rows-removed_rows.size()];
    size_t new_id = 0;
    for (size_t i=0; i<n_org_rows; i++) {
        if (find(removed_rows.begin(), removed_rows.end(), i)!=removed_rows.end()) continue;
        (*eqsRHS)[new_id] = org_eqsRHS[i];
        (*eqsX)[new_id] = org_eqsX[i];
        map_solved_orgEqs[new_id] = i;
        new_id++;
    }
    run.stop();
    std::cout << run.elapsed() << " s" << std::endl;
}

void mapSolvedXToOriginalX(double *eqsX, size_t dim, map<INDEX_TYPE,INDEX_TYPE> &map_solved_orgEqs, double *org_eqsX)
{
    for (size_t i=0; i<dim; i++) {
    	INDEX_TYPE idx (map_solved_orgEqs[i]);
    	double tmp = eqsX[i];
        org_eqsX[idx] = tmp;
    }
}
#endif

//-----------------------------------------------------------------------------
// Test 1
//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
//output current setting
    std::cout << "##### CURRENT SETTING #####" << std::endl;
#ifdef LIS
    std::cout << "- Linear Solver: LIS" << std::endl;
#else
    std::cout << "- Linear Solver: Tom" << std::endl;
#endif
#ifdef USE_EIGEN
    std::cout << "- Matrix Library: Eigen" << std::endl;
#else
    std::cout << "- Matrix Library: CRSMatrix" << std::endl;
#endif
    std::cout << "- Option: ";
#ifdef NO_SOLVER
    std::cout << "Skip solver / ";
#endif
#ifdef OUTPUT_VTK
    std::cout << "Output VTK / ";
#endif
    std::cout << std::endl;
    std::cout << "###########################" << std::endl;

#ifdef LIS
    lis_initialize(&argc, &argv);
#endif
    double eps (1.0e-10);
    unsigned steps (10000);
    if (argc > 2)
        sscanf(argv[2], "%lf", &eps);
    if (argc > 3)
        sscanf(argv[3], "%d", &steps);
#ifdef _OPENMP
    int nthreads = 1;
    if (argc > 4)
        sscanf(argv[4], "%d", &nthreads);
    omp_set_num_threads (nthreads);
#endif
#ifdef LIS
    MathLib::LIS_option option;
    option.ls_method = 1;
    option.ls_precond = 0;
    option.ls_extra_arg = "";
    option.ls_max_iterations = steps;
    option.ls_error_tolerance = eps;
    if (argc > 5)
        sscanf(argv[5], "%d", &option.ls_method);
    if (argc > 6)
        sscanf(argv[6], "%d", &option.ls_precond);
#endif
#ifdef _OPENMP
    cout << "->Start OpenMP parallelization with " << omp_get_max_threads() << " threads" << endl;
#if 1
    {
      int a[1000000];
      RunTimeTimer run_timer;
      CPUTimeTimer cpu_timer;
      run_timer.start();
      cpu_timer.start();
      double time0 = omp_get_wtime();
      #pragma omp parallel for
      for (int i=0; i<1000000; i++)
        a[i] = 0;
      run_timer.stop();
      cpu_timer.stop();
      double time1 = omp_get_wtime();
      cout.setf(std::ios::scientific,std::ios::floatfield);
      cout.precision(12);
      cout << "CPU time = " << cpu_timer.elapsed() << endl;
      cout << "Run time = " << run_timer.elapsed() << endl;
      //cout << "Wall time= " << (time1-time0) << endl;
    }
#endif
#endif
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
    msh->construct();

    //set material properties
    const double K = 1.e-8;
    const double S = .0;

    //set Dirichlet BC
    vector<IndexValue> list_dirichlet_bc;
    setDirichletBC_Case1(msh, list_dirichlet_bc);

    //-- construct EQS -----------------------------------------------
    //prepare EQS
    const size_t dim_eqs = msh->getNumberOfNodes();
    MathLib::SparseTableCRS<INDEX_TYPE>* crs = FemLib::generateSparseTableCRS<INDEX_TYPE>(msh);
    //FemLib::outputSparseTableCRS(crs);
#ifdef USE_EIGEN
    Eigen::MappedSparseMatrix<double, Eigen::RowMajor> eqsA(dim_eqs, dim_eqs, crs->nonzero, crs->row_ptr, crs->col_idx, crs->data);
#else

#ifdef _OPENMP
#ifdef LIS
    MathLib::CRSMatrixOpenMP<double,int> eqsA(static_cast<unsigned>(crs->dimension), crs->row_ptr, crs->col_idx, crs->data, nthreads);
#else
    MathLib::CRSMatrixOpenMP<double,unsigned> eqsA(static_cast<unsigned>(crs->dimension), crs->row_ptr, crs->col_idx, crs->data, nthreads);
#endif
#else
    MathLib::CRSMatrix<double, INDEX_TYPE> eqsA(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
//    MathLib::CRSMatrixDiagPrecond eqsA(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
#endif

#endif
    double* eqsX(new double[dim_eqs]);
    double* eqsRHS(new double[dim_eqs]);
    for (size_t i=0; i<dim_eqs; i++) eqsX[i] = 0.0;
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
    RunTimeTimer run_timer_assembly;
    CPUTimeTimer cpu_timer_assembly;
    run_timer_assembly.start();
    cpu_timer_assembly.start();
    double time0 = omp_get_wtime();

    #ifdef _OPENMP
    #pragma omp parallel for default(none) private(ele, dof_map, pt, nodes_x, nodes_y, nodes_z, a, b, c, local_K) shared(msh, eqsA)
    #endif
    for (long i_ele=0; i_ele<static_cast<long>(n_ele); i_ele++) {
        ele = static_cast<MeshLib::Triangle*>(msh->getElemenet(i_ele));
        // setup element information
        // dof
        for (size_t i=0; i<ele->getNumberOfNodes(); i++)
            dof_map[i] = ele->getNodeID(i);
        // xyz
        for (size_t i=0; i<n_ele_nodes; i++) {
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
        for (size_t i=0; i<n_ele_nodes; i++)
            for (size_t j=0; j<i; j++)
                local_K(i,j) = local_K(j,i);

        // add into global EQS
        for (size_t i=0; i<n_ele_nodes; i++) {
            for (size_t j=0; j<n_ele_nodes; j++) {
#ifdef USE_EIGEN
              #pragma omp atomic
              eqsA.coeffRef(dof_map[i], dof_map[j]) += local_K(i,j);
#else
              eqsA.addValue(dof_map[i], dof_map[j], local_K(i,j));
#endif
            }
        }
    }
    double time1 = omp_get_wtime();
    run_timer_assembly.stop();
    cpu_timer_assembly.stop();

    //MathLib::EigenTools::outputEQS("eqs1.txt", eqsA, eqsX, eqsRHS);

    cout << "->apply BC" << endl;
    CPUTimeTimer cpu_timer3;
    RunTimeTimer run_timer3;
    run_timer3.start();
    cpu_timer3.start();

    //apply Dirichlet BC
#if 1
#ifdef USE_EIGEN
    for (size_t i=0; i<list_dirichlet_bc.size(); i++) {
        IndexValue &bc = list_dirichlet_bc.at(i);
        const size_t id = bc.id;
        const double val = bc.val;
        MathLib::EigenTools::setKnownXi(eqsA, eqsRHS, id, val);
    }
#else
    double *org_eqsX = eqsX;
    double *org_eqsRHS = eqsRHS;
    map<INDEX_TYPE,INDEX_TYPE> map_solved_orgEqs;
    setKnownXi_ReduceSizeOfEQS(list_dirichlet_bc, eqsA, org_eqsRHS, org_eqsX, &eqsRHS, &eqsX, map_solved_orgEqs);
#endif
    run_timer3.stop();
    cpu_timer3.stop();
#endif

    //apply ST
    //MathLib::EigenTools::outputEQS("eqs2.txt", eqsA, eqsX, eqsRHS);

#ifndef USE_EIGEN
    // output matrix
//    std::cout << "writing matrix to matrix.bin ... " << std::flush;
//    std::ofstream out ("matrix.bin", std::ios::binary);
//    CS_write (out, eqsA.getNRows(), eqsA.getRowPtrArray(), eqsA.getColIdxArray(), eqsA.getEntryArray());
//    std::cout << "ok" << std::endl;
#endif

//	std::string fname_fem_out (strMeshFile);
//	fname_fem_out = fname_fem_out.substr(0,fname_fem_out.length()-4);
//	fname_fem_out += "_fem.bin";
//	std::ofstream os (fname_fem_out.c_str(), std::ios::binary);
//	if (os) {
//		std::cout << "writing FEM matrix to " << fname_fem_out << " ... " << std::flush;
//		CS_write(os, eqsA.getNRows(), eqsA.getRowPtrArray(), eqsA.getColIdxArray(), eqsA.getEntryArray());
//		std::cout << "ok" << std::endl;
//		os.close();
//		std::ofstream out_rhs("rhs.dat");
//		if (out_rhs) {
//			const unsigned n_entries(eqsA.getNRows());
//			for (unsigned k(0); k<n_entries; k++) {
//				out_rhs << eqsRHS[k] << std::endl;
//			}
//			out_rhs.close();
//		}
//	} else {
//		std::cout << "could not open file " << fname_fem_out << " for writing" << std::endl;
//	}

    //-- solve EQS -----------------------------------------------
    //set up
    cout << "->solve EQS" << endl;
    CPUTimeTimer cpu_timer2;
    RunTimeTimer run_timer2;
    run_timer2.start();
    cpu_timer2.start();

    //omp_set_num_threads (1);

#ifndef NO_SOLVER
#ifdef LIS

#ifdef USE_EIGEN
    MathLib::solveWithLis(crs, eqsX, eqsRHS, option);
#else
    crs = new MathLib::SparseTableCRS<int>();
    crs->nonzero = eqsA.getNNZ();
    crs->dimension = eqsA.getNRows();
    crs->row_ptr = (int*)eqsA.getRowPtrArray();
    crs->col_idx = (int*)eqsA.getColIdxArray();
    crs->data = (double*)eqsA.getEntryArray();
    //FemLib::outputSparseTableCRS(crs);

    MathLib::solveWithLis(crs, eqsX, eqsRHS, option);

    mapSolvedXToOriginalX(eqsX, crs->dimension, map_solved_orgEqs, org_eqsX);
    double *temp_x = eqsX;
    eqsX = org_eqsX;
    org_eqsX = temp_x;
#endif
#else
#ifndef USE_EIGEN
//    eqsA.calcPrecond();
    std::cout << "solving system of " << eqsA.getNRows() << " linear equations (nnz = " << eqsA.getNNZ() << ") " << std::flush;
#ifndef _OPENMP
    std::cout << " with sequential solver" << std::endl;
    MathLib::CG(&eqsA, eqsRHS, eqsX, eps, steps);
    std::cout << "MathLib::CG converged within " << steps << ", residuum is " << eps << std::endl;
#else
    if (nthreads == 1) {
    	std::cout << " with sequential solver" << std::endl;
    	MathLib::CG(&eqsA, eqsRHS, eqsX, eps, steps);
    } else {
    	std::cout << " with OpenMP parallelized solver" << std::endl;
    	MathLib::CGParallel (&eqsA, eqsRHS, eqsX, eps, steps);
    }
	std::cout << "MathLib::CGParallel converged within " << steps << ", residuum is " << eps << std::endl;
#endif
    mapSolvedXToOriginalX(eqsX, eqsA.getNRows(), map_solved_orgEqs, org_eqsX);
	double *temp_x = eqsX;
	eqsX = org_eqsX;
	org_eqsX = temp_x;
#endif
#endif
#endif
    run_timer2.stop();
    cpu_timer2.stop();

    run_timer.stop();
    cpu_timer.stop();
    cout.setf(std::ios::scientific,std::ios::floatfield);
    cout.precision(12);
    cout << "== Simulation time ==" << endl;
    cout << "Total simulation:" << endl;
    cout << "CPU time = " << cpu_timer.elapsed() << endl;
    cout << "Run time = " << run_timer.elapsed() << endl;
    cout << "---------------------" << endl;
    cout << "Assembly:" << endl;
    cout << "CPU time = " << cpu_timer_assembly.elapsed() << endl;
    cout << "Run time = " << run_timer_assembly.elapsed() << endl;
//    cout << "Wall time= " << (time1-time0) << endl;
    cout << "---------------------" << endl;
    cout << "Apply BC:" << endl;
    cout << "CPU time = " << cpu_timer3.elapsed() << endl;
    cout << "Run time = " << run_timer3.elapsed() << endl;
    cout << "---------------------" << endl;
    cout << "Linear solver:" << endl;
    cout << "CPU time = " << cpu_timer2.elapsed() << endl;
    cout << "Run time = " << run_timer2.elapsed() << endl;

    // output results
#ifdef OUTPUT_VTK
    cout << "->output results" << endl;
    std::vector<MeshLib::NodalScalarValue> nodalValues;
    string str = "Head";
    MeshLib::NodalScalarValue temp("Head", eqsX);
    nodalValues.push_back(temp);
    MeshLib::MeshIOLegacyVtk4Simulation::WriteAsciiFile("output.vtk", *msh, 1, 1.0, nodalValues);
#endif
    //release memory
#ifdef LIS
    lis_finalize();
#endif
    destroyStdVectorWithPointers(vec_mesh);
    //delete [] crs->row_ptr;
    //delete [] crs->col_idx;
    ////delete [] crs->data;
#ifdef USE_EIGEN
    delete crs;
#else
    delete[] org_eqsX;
    delete[] org_eqsRHS;
#endif
    delete [] eqsX;
    delete [] eqsRHS;

    return 0;
}
