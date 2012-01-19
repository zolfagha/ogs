/*
 * FEMTest-OGS-6-CRS-OGS-6-CG-Perm.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: TF
 */


#include <iostream>
#include <limits>
#include <algorithm>
#include <map>

#include <Eigen>

//#include "LinAlg/Sparse/CRSMatrixOpenMP.h"
//#include "LinAlg/Sparse/CRSMatrixDiagPrecond.h"
//#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Sparse/SparseTableCRS.h"

#include "sparse.h"

#include "LinAlg/Solvers/CG.h"
#include "Mesh.h"
#include "IO/MeshIOOGSAscii.h"
#include "IO/MeshIOLegacyVtk.h"
#include "Tools/MeshGenerator.h"
#include "CPUTimeTimer.h"
#include "RunTimeTimer.h"
#include "MeshSparseTable.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/AdjMat.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReordered.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/Cluster.h"
#include "LinAlg/Sparse/CRSMatrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define INDEX_TYPE unsigned

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
    if (abs(pt[0]-x_min)<std::numeric_limits<double>::epsilon()) {
      list_dirichlet_bc.push_back(IndexValue(i, head_left));
    } else if (abs(pt[0]-x_max)<std::numeric_limits<double>::epsilon()) {
      list_dirichlet_bc.push_back(IndexValue(i, head_right));
    }
  }
}

#ifndef USE_EIGEN
void setKnownXi_ReduceSizeOfEQS(std::vector<IndexValue> &list_dirichlet_bc, MathLib::CRSMatrix<
				double, INDEX_TYPE> &eqsA, double* org_eqsRHS, double* org_eqsX, double** eqsRHS,
				double** eqsX, std::map<INDEX_TYPE, INDEX_TYPE> &map_solved_orgEqs)
{
	const size_t n_org_rows = eqsA.getNRows();
	std::vector<INDEX_TYPE> removed_rows(list_dirichlet_bc.size());
	std::cout << "\t[BC] (transpose matrix) ... " << std::flush;
	RunTimeTimer run_trans;
	run_trans.start();
	MathLib::CRSMatrix<double, INDEX_TYPE>* transpose_mat (eqsA.getTranspose());
	run_trans.stop();
	std::cout << run_trans.elapsed() << " s" << std::endl;

	INDEX_TYPE const*const row_ptr (transpose_mat->getRowPtrArray());
	INDEX_TYPE const*const col_idx (transpose_mat->getColIdxArray());
	double const*const data (transpose_mat->getEntryArray());

	std::cout << "\t[BC] modifying rhs and solution vector ... " << std::flush;
	RunTimeTimer run;
	run.start();
	for (size_t i = 0; i < list_dirichlet_bc.size(); i++) {
		IndexValue &bc = list_dirichlet_bc.at(i);
		const size_t id = bc.id;
		const double val = bc.val;
		removed_rows.at(i) = id;

		//b_i -= A(i,k)*val, i!=k
//		for (size_t j = 0; j < eqsA.getNCols(); j++)
//			org_eqsRHS[j] -= eqsA.getValue(j, id) * val;

		const INDEX_TYPE end(row_ptr[id+1]);
		for (INDEX_TYPE k(row_ptr[id]); k<end; k++) {
			const INDEX_TYPE j(col_idx[k]);
			org_eqsRHS[j] -= data[k] * val;
		}

		//b_k = A_kk*val
		org_eqsRHS[id] = val; //=eqsA(id, id)*val;
		org_eqsX[id] = val; //=eqsA(id, id)*val;
	}
	run.stop();
	std::cout << run.elapsed() << " s" << std::endl;
	delete transpose_mat;

	std::cout << "\t[BC] erasing rows and columns from matrix ... " << std::flush;
	run.start();
	//remove rows and columns
	eqsA.eraseEntries(removed_rows.size(), &removed_rows[0]);
	run.stop();
	std::cout << run.elapsed() << " s" << std::endl;

	//remove X,RHS
	std::cout << "\t[BC] create mapping ... " << std::flush;
	run.start();
	(*eqsX) = new double[n_org_rows - removed_rows.size()];
    (*eqsRHS) = new double[n_org_rows-removed_rows.size()];
    size_t new_id = 0;
    for (size_t i=0; i<n_org_rows; i++) {
        if (std::find(removed_rows.begin(), removed_rows.end(), static_cast<unsigned>(i))!=removed_rows.end()) continue;
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
        org_eqsX[map_solved_orgEqs[i]] = eqsX[i];
    }
}
#endif

int main(int argc, char *argv[])
{
    double eps (1.0e-6);
    unsigned steps (10000);
    if (argc > 2)
        sscanf(argv[2], "%lf", &eps);
    if (argc > 3)
        sscanf(argv[3], "%d", &steps);

#ifdef _OPENMP
	omp_set_num_threads(1);
#endif
//    std::cout << "->Start OpenMP parallelization with " << omp_get_max_threads() << " threads" << std::endl;

    // *** setup a problem
    // set mesh
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

    MeshLib::UnstructuredMesh *msh = static_cast<MeshLib::UnstructuredMesh*>(vec_mesh.at(0));
    if (msh->getElemenet(0)->getElementType()!=MeshLib::ElementType::TRIANGLE) {
        std::cout << "Only triangle elements are supported! " << endl;
        return 0;
    }
    msh->construct();

    // set material properties
//    const double K = 1.e-8;
//    const double S = .0;

    //set Dirichlet BC
    std::vector<IndexValue> list_dirichlet_bc;
    setDirichletBC_Case1(msh, list_dirichlet_bc);

    // *** construct EQS
    // prepare EQS
    const size_t dim_eqs = msh->getNumberOfNodes();
    MathLib::SparseTableCRS<INDEX_TYPE>* crs = FemLib::generateSparseTableCRS<INDEX_TYPE>(msh);

    MathLib::CRSMatrixReordered eqsA(static_cast<unsigned>(crs->dimension), crs->row_ptr, crs->col_idx, crs->data);

    double* eqsX(new double[dim_eqs]);
    double* eqsRHS(new double[dim_eqs]);
    for (size_t i=0; i<dim_eqs; i++) eqsX[i] = 0.0;
    for (size_t i=0; i<dim_eqs; i++) eqsRHS[i] = 0.0;

    // assembly EQS
    const size_t n_ele = msh->getNumberOfElements();
    Eigen::Matrix3d local_K = Eigen::Matrix3d::Zero();
    const size_t n_ele_nodes = 3;
    double nodes_x[n_ele_nodes], nodes_y[n_ele_nodes], nodes_z[n_ele_nodes];
    double a[n_ele_nodes];
    double b[n_ele_nodes];
    double c[n_ele_nodes];
    size_t dof_map[n_ele_nodes];
    MeshLib::Triangle *ele;
    double pt[3];

    std::cout << "-> assembly EQS" << std::endl;

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
              eqsA.addValue(dof_map[i], dof_map[j], local_K(i,j));
            }
        }
    }

    cout << "-> apply BC" << endl;
    CPUTimeTimer cpu_timer3;
    RunTimeTimer run_timer3;
    run_timer3.start();
    cpu_timer3.start();
    // apply Dirichlet BC
    double *org_eqsX = eqsX;
    double *org_eqsRHS = eqsRHS;
    map<INDEX_TYPE,INDEX_TYPE> map_solved_orgEqs;
    setKnownXi_ReduceSizeOfEQS(list_dirichlet_bc, eqsA, org_eqsRHS, org_eqsX, &eqsRHS, &eqsX, map_solved_orgEqs);
    run_timer3.stop();
    cpu_timer3.stop();
    std::cout << "\t[BC] sum CPU time = " << run_timer3.elapsed() << std::endl;
    std::cout << "\t[BC] sum Run time = " << cpu_timer3.elapsed() << std::endl;
    // writing system of linear equations to file for external solver
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

    bool verbose(true);
    // *** permute matrix and rhs with nested dissection approach
    // create time measurement objects
	RunTimeTimer run_timer_nd;
    // calculate the nested dissection reordering
	if (verbose) {
		std::cout << "-> calculating nested dissection reordering of matrix ... " << std::flush;
	}
	run_timer_nd.start();
	const size_t n(eqsA.getNRows());
	MathLib::Cluster cluster_tree(n, const_cast<unsigned*>(eqsA.getRowPtrArray()),
					const_cast<unsigned*>(eqsA.getColIdxArray()));
	unsigned *op_perm(new unsigned[n]);
	unsigned *po_perm(new unsigned[n]);
	for (unsigned k(0); k < n; k++)
		op_perm[k] = po_perm[k] = k;
	cluster_tree.createClusterTree(op_perm, po_perm, 1000);
	run_timer_nd.stop();
	if (verbose) {
		std::cout << run_timer_nd.elapsed() << " s" << std::endl;
	}
	// applying the nested dissection reordering to matrix
	RunTimeTimer run_timer_apl_nd;
	if (verbose) {
		std::cout << "-> applying nested dissection reordering to FEM matrix ... " << std::flush;
	}
	run_timer_apl_nd.start();
	eqsA.reorderMatrix(op_perm, po_perm);
	run_timer_apl_nd.stop();
	if (verbose)
		std::cout << run_timer_apl_nd.elapsed() << " s" << std::endl;

	// applying the nested dissection reordering to rhs
	RunTimeTimer run_timer_apl_nd_rhs;
	run_timer_apl_nd_rhs.start();
	if (verbose) {
		std::cout << "-> applying nested dissection reordering to rhs ... " << std::flush;
	}
	double *tmp_rhs(new double[eqsA.getNRows()]);
	for (size_t k(0); k<n; k++) tmp_rhs[k] = eqsRHS[op_perm[k]];
	for (size_t k(0); k<n; k++) eqsRHS[k] = tmp_rhs[k];
	delete [] tmp_rhs;
	run_timer_apl_nd_rhs.stop();
	if (verbose) {
		std::cout << run_timer_apl_nd_rhs.elapsed() << " s" << std::endl;
	}

    // *** solve EQS
    std::cout << "-> solve EQS" << std::endl;
    CPUTimeTimer cpu_timer2;
    RunTimeTimer run_timer2;
    run_timer2.start();
    cpu_timer2.start();
//    eqsA.calcPrecond();
    std::cout << "-> solving system of " << eqsA.getNRows() << " linear equations (nnz = " << eqsA.getNNZ() << ") " << std::flush;
    std::cout << " with sequential solver" << std::endl;
    MathLib::CG(&eqsA, eqsRHS, eqsX, eps, steps);
    std::cout << "MathLib::CG converged within " << steps << ", residuum is " << eps << std::endl;

    // applying the nested dissection reordering to the solution
	if (verbose) {
		std::cout << "-> reverse the nested dissection reordering to solution x ... " << std::flush;
	}
    double *tmp_x(new double[eqsA.getNRows()]);
	for (size_t k(0); k<n; k++) tmp_x[k] = eqsX[po_perm[k]];
	for (size_t k(0); k<n; k++) eqsX[k] = tmp_x[k];
	delete [] tmp_x;
	if (verbose) {
		std::cout << "done" << std::endl;
	}

//    mapSolvedXToOriginalX(eqsX, crs->dimension, map_solved_orgEqs, org_eqsX);
    mapSolvedXToOriginalX(eqsX, eqsA.getNRows(), map_solved_orgEqs, org_eqsX);

//	std::ofstream vec_out ("vector_perm.txt");
//	if (vec_out) {
//		vec_out.precision(3);
//		for (size_t k(0); k<n; k++)
//			vec_out << eqsX[k] << std::endl;
//	}

	double *temp_x = eqsX;
	eqsX = org_eqsX;
	org_eqsX = temp_x;
	run_timer2.stop();
    cpu_timer2.stop();

    run_timer.stop();
    cpu_timer.stop();
    std::cout.setf(std::ios::scientific,std::ios::floatfield);
    std::cout.precision(12);
    std::cout << "== Simulation time ==" << std::endl;
    std::cout << "Total simulation:" << std::endl;
    std::cout << "CPU time = " << run_timer.elapsed() << std::endl;
    std::cout << "Run time = " << cpu_timer.elapsed() << std::endl;
    std::cout << "Apply BC:" << std::endl;
    std::cout << "CPU time = " << run_timer3.elapsed() << std::endl;
    std::cout << "Run time = " << cpu_timer3.elapsed() << std::endl;
    std::cout << "NestedDissection: " << std::endl;
    std::cout << "Run time = " << run_timer_nd.elapsed() + run_timer_apl_nd.elapsed() + run_timer_apl_nd_rhs.elapsed() << std::endl;
    std::cout << "ND + lin solver: " << std::endl;
    std::cout << "CPU time = " << run_timer2.elapsed() + run_timer_nd.elapsed() + run_timer_apl_nd.elapsed() + run_timer_apl_nd_rhs.elapsed() << std::endl;
    std::cout << "Linear solver:" << std::endl;
    std::cout << "CPU time = " << run_timer2.elapsed() << std::endl;
    std::cout << "Run time = " << cpu_timer2.elapsed() << std::endl;

    // output results
    cout << "-> output results" << endl;
    std::vector<MeshLib::NodalScalarValue> nodalValues;
    string str = "Head";
    MeshLib::NodalScalarValue temp("Head", eqsX);
    nodalValues.push_back(temp);
    MeshLib::MeshIOLegacyVtk4Simulation::WriteAsciiFile("output.vtk", *msh, 1, 1.0, nodalValues);

    //release memory
    destroyStdVectorWithPointers(vec_mesh);
    delete[] org_eqsX;
    delete[] org_eqsRHS;
    delete [] eqsX;
    delete [] eqsRHS;

    return 0;
}

