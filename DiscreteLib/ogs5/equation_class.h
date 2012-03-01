/**************************************************************************
   Task: Sparse matrix and linear equation solver
      Design and program by WW
   Programing:
   10/2007 WW/
**************************************************************************/
#ifndef eqs_class_INC
#define eqs_class_INC

#include <cmath>
#include <iostream>
#include <vector>

#include "matrix_class.h"
#include "DenseMatrix.h"

namespace OGS5
{

class CNumerics;
class CPARDomain;



class Linear_EQS
{
public:
	Linear_EQS(const SparseTable &sparse_table, const long dof, bool messg = true);
#ifdef USE_MPI
	Linear_EQS(const long size);
#endif
	~Linear_EQS();

	void ConfigNumerics( CNumerics* m_num, const long n = 0);

    //
    void Initialize();
    void Clean();

#if defined(USE_MPI)
	int Solver(double* xg, const long n);
	double GetCPUtime() const { return cpu_time;  }
#else
	int Solver();
#endif

	//
	// Access to the members
	void SetDOF(const int dof_n)
	{
		A->SetDOF(dof_n);
	}
	void SetKnownX_i(const long i, const double x_i);
	double X(const long i) const {return x[i]; }
    double* getX() {return x;};
	double RHS(const long i) const {return b[i]; }
    double* getRHS() {return b;};
	double NormX();
	double NormRHS() { return bNorm; }
#if defined(USE_MPI)
	int DOF() { return A->Dof(); }
	long Size() { return A->Size(); }
	void SetDomain(CPARDomain* a_dom) {dom = a_dom; }
#endif
	// Write
	void Write(std::ostream &os = std::cout);
	void WriteRHS(std::ostream &os = std::cout);
	void WriteX(std::ostream &os = std::cout);
	void Write_BIN(std::ostream &os);

    CSparseMatrix* getA() {return A;};

private:                                          // Dot not remove this!
	CSparseMatrix* A;
	double* b;
	double* x;
	double* prec_M;
	//
#if defined(USE_MPI)
	CPARDomain* dom;
	// WW
	double* border_buffer0;
	double* border_buffer1;
	double cpu_time;
	//
	double dot (const double* xx,  const double* yy, const long n);
	inline void MatrixMulitVec(double* xx,  double* yy);
	inline void TransMatrixMulitVec(double* xx,  double* yy);
#endif
	//
	std::string solver_name;
	std::string precond_name;
	// Buffer
	std::vector<double*> f_buffer;
	// Controls
	int precond_type;
	int solver_type;
	bool message;
	int iter, max_iter;
	double tol, bNorm, error;
	long size_global;
	long size_A;
	// Operators
	double dot (const double* xx,  const double* yy);
	inline double Norm(const double* xx)  { return sqrt(dot(xx, xx)); }
	inline bool CheckNormRHS(const double normb_new);
	/// GMRES. 30.06.2010. WW
	/// GMRES H matrix
	mutable Matrix H;
	int m_gmres;                          /// number of columns of H matrix
	void Update(double* x, int k, Matrix &h, double* s);
	void Get_Plane_Rotation(double &dx, double &dy, double &cs, double &sn);
	void Set_Plane_Rotation(double &dx, double &dy, double &cs, double &sn);
	//
	void Message();

    // Preconditioner;
    void Precond(double* vec_s, double* vec_r);
    void TransPrecond(double* vec_s, double* vec_r);
    //#if defined(USE_MPI)
    void Precond_Jacobi(const double* vec_s, double* vec_r);
    //#endif
#if defined(USE_MPI)
    int CG(double* xg, const long n);
    int BiCG(double* xg, const long n);   //02.2010. WW
    int BiCGStab(double* xg, const long n);
    int CGS(double* xg, const long n);
#else
    int CG();
    int BiCG();                           //02.2010. WW
    int BiCGStab();
    int Gauss() {return -1; }
    int QMRCGStab() {return -1; }
    int CGNR() {return -1; }
    int CGS();
    int Richardson() {return -1; }
    int JOR() {return -1; }
    int SOR() {return -1; }
    int AMG1R5() {return -1; }
    int UMF() {return -1; }
    int GMRES();
#endif
    //
    void ComputePreconditioner();
    void ComputePreconditioner_Jacobi();
    void ComputePreconditioner_ILU() {       }
};
}

#endif
