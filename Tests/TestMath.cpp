
#include <gtest/gtest.h>

#include <valarray>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"
#include "MathLib/Nonlinear/Picard.h"

#include "TestUtil.h"

using namespace MathLib;



/**
 * \brief Differentiable function
 */
template <class T_D0, class T_D1>
class IDifferentiableFunction
{
public:
    virtual ~IDifferentiableFunction() {};
    virtual void eval0(T_D0 &x, T_D0 &r) = 0;
    virtual void eval1(T_D0 &x, T_D1 &j) = 0;
};


// Example problems ------------------------------------------------------------------
// f(x) = x*x -4 = 0
// x = 2,-2
class NL1_NR : public IDifferentiableFunction<double,double>
{
public:
	// r = f(x)
	void eval0(double &x, double &r) { r = x*x-4.; }
	// j = f'(x)
	void eval1(double &x, double &j) { j = 2*x; }
};

class NL1_PC
{
public:
	// x = g(x) = x - 1/f'(x) * f(x)
	void eval0(double &x, double &x_new) { x_new = (x*x+4.)/(2.*x); }
};


// x*y - y = 0
// y*y - 1 = 0
// x=1, y=1,-1
class NL2_NR : public IDifferentiableFunction<std::valarray<double>,MathLib::Matrix<double> >
{
public:
	void eval0(std::valarray<double> &x, std::valarray<double> &r)
	{
		r[0] = x[0]*x[1]-x[1];
		r[1] = x[1]*x[1]-1.0;
	}
	void eval1(std::valarray<double> &x, MathLib::Matrix<double> &j)
	{
		j(0,0) = x[1];
		j(0,1) = x[0]-1.0;
		j(1,0) = .0;
		j(1,1) = 2.*x[1];
	}
};

template<class T_DIFF_FUNCTION>
class NL2_PC
{
	MathLib::DenseLinearEquations solver;
	std::valarray<double> b;
	std::valarray<double> tmp;
	T_DIFF_FUNCTION nr;
public:
	NL2_PC(size_t n)
	{
		solver.create(n);
		b.resize(n);
		tmp.resize(n);
	}
	// x = g(x) = x - f'(x)^-1 * f(x)
	// f'*g = f'*x - f
	// (if f=A*x+b, A*g=b)
	void eval0(std::valarray<double> &x, std::valarray<double> &x_new)
	{
		MathLib::Matrix<double> *j = solver.getA();
		solver.reset();
		b = 0;
		tmp = 0;

		nr.eval1(x, *j);
		nr.eval0(x, b);

		b*=-1;
		j->axpy(1.0, &x[0], .0, &tmp[0]);
		b+=tmp;
		for (size_t i=0; i<x.size(); i++)
			solver.setRHS(i, b[i]);

		solver.solve();

    	double *u = solver.getX();
		for (size_t i=0; i<x.size(); i++)
			x_new[i] = u[i];
	}

};

class NL2_PC2
{
	MathLib::DenseLinearEquations solver;
	std::valarray<double> b;
public:
	NL2_PC2(size_t n)
	{
		solver.create(n);
		b.resize(n);
	}
	// if f=A*x+b, A*g=b
	void eval0(std::valarray<double> &x, std::valarray<double> &x_new)
	{
		MathLib::Matrix<double> *j = solver.getA();
		solver.reset();
		b = 0;

		(*j)(0,0) = 3;
		(*j)(0,1) = -1.0;
		(*j)(1,0) = 2*x[0];
		(*j)(1,1) = -1;
		b[0] = -2;
		b[1] = 0.;
        //(*j)(0,0) = x[1];
        //(*j)(0,1) = -1.0;
        //(*j)(1,0) = .0;
        //(*j)(1,1) = x[1];
        //b[0] = .0;
        //b[1] = 1.;
		for (size_t i=0; i<x.size(); i++)
			solver.setRHS(i, b[i]);

		solver.solve();

    	double *u = solver.getX();
		for (size_t i=0; i<x.size(); i++)
			x_new[i] = u[i];
	}

};


class NL3_NR : public IDifferentiableFunction<std::valarray<double>, MathLib::CRSMatrix<double, signed> >
{
public:
  NL3_NR()
  {
	  const size_t n = 2;
	  MathLib::RowMajorSparsity sparse(n);
	  sparse[0].insert(0);
	  sparse[0].insert(1);
	  sparse[1].insert(0);
	  sparse[1].insert(1);
	  _linear_solver.create(n, &sparse);
  }

  void eval0(std::valarray<double> &x, std::valarray<double> &r)
  {
	  r[0] = x[0]*x[1]-x[1];
	  r[1] = x[1]*x[1]-1.0;
  }

  void eval1(std::valarray<double> &x, MathLib::CRSMatrix<double, signed> &j)
  {
	  j.setValue(0, 0, x[1]);
	  j.setValue(0, 1, x[0]-1.0);
	  j.setValue(1, 0, .0);
	  j.setValue(1, 1,  2.*x[1]);
  }

  MathLib::CRSLisSolver* getLinearSolver() {return &_linear_solver;};


private:
  MathLib::CRSLisSolver _linear_solver;
};

//------------------------------------------------------------------

TEST(Math, NonlinearNR_double)
{
	NL1_NR f1;
	double x0 = 6.0;
	double x = .0;
	NewtonRaphsonMethod nr;
	nr.solve(f1, x0, x);
	ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearNR_dense)
{
	typedef std::valarray<double> MyVector;

	NL2_NR f2;
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	NewtonRaphsonMethod nr;
	nr.solve(f2, x0, x);

	double my_expect[] = {1., 1.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearNR_sparse)
{
	typedef std::valarray<double> MyVector;
	typedef MathLib::CRSMatrix<double, signed> MyMatrix;
	typedef NRCheckConvergence<MyVector,NRErrorNorm1DX > MyConverge;
	typedef DxSolverVector<MyVector,MyMatrix> MyDxSolver;

	NL3_NR f;
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	MyVector r(2), dx(2);
	MyDxSolver dx_solver(f.getLinearSolver());
	MyMatrix* j = f.getLinearSolver()->getA();

	NewtonRaphsonMethod nr;
	nr.solve<NL3_NR,MyVector,MyMatrix,MyDxSolver,MyConverge>(f, x0, x, r, *j, dx, dx_solver);

	double my_expect[] = {1., 1.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearPicard_double)
{
	NL1_PC f1;
	double x0 = 6.0;
	double x = .0;
	PicardMethod nr;
	nr.solve(f1, x0, x);
	ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearPicard_dense1)
{
	typedef std::valarray<double> MyVector;

	NL2_PC<NL2_NR> f2(2);
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	PicardMethod nr;
	nr.solve(f2, x0, x);

	double my_expect[] = {1., 1.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearPicard_dense2)
{
	typedef std::valarray<double> MyVector;

	NL2_PC2 f2(2);
	MyVector x0(1.1, 2);
	MyVector x(0.0, 2);
	PicardMethod nr;
	nr.solve(f2, x0, x);

	double my_expect[] = {-0.5, 0.5};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
//    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}


TEST(Math, Matrix_transposeAndMultiply)
{
    MathLib::Matrix<double> matA(1,2);
    MathLib::Matrix<double> matB(1,2);
    MathLib::Matrix<double> matC(2,2);

    matA(0,0) = 1.0;
    matA(0,1) = 2.0;
    matB(0,0) = 3.0;
    matB(0,1) = 4.0;
    matC = .0;
    matA.transposeAndMultiply(matB, matC);

    ASSERT_EQ(matC(0,0), 3.0);
    ASSERT_EQ(matC(0,1), 4.0);
    ASSERT_EQ(matC(1,0), 6.0);
    ASSERT_EQ(matC(1,1), 8.0);
}


TEST(Math, MatrixAddSubMatrix)
{
    Matrix<double> _m(4,4);
    _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0; 
    _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
    _m(2,2) = 4.0; _m(2,3) = -1.0;
    _m(3,3) = 4.0;
    for (size_t i=0; i<4; i++)
        for (size_t j=0; j<i; j++) _m(i,j) = _m(j,i);

    DenseLinearEquations eqs;
    eqs.create(10);

    for (size_t i=0; i<3; i++) {
        std::vector<long> pos(4);
        for (size_t j=0; j<4; j++)
            pos[j] = i*3 + j;
        eqs.addAsub(pos, _m);
    }

    double expected_m[] = {
         4, -1, -2, -1,  0,  0,  0,  0,  0,  0,
        -1,  4, -1, -2,  0,  0,  0,  0,  0,  0,
        -2, -1,  4, -1,  0,  0,  0,  0,  0,  0,
        -1, -2, -1,  8, -1, -2, -1,  0,  0,  0,
         0,  0,  0, -1,  4, -1, -2,  0,  0,  0,
         0,  0,  0, -2, -1,  4, -1,  0,  0,  0,
         0,  0,  0, -1, -2, -1,  8, -1, -2, -1,
         0,  0,  0,  0,  0,  0, -1,  4, -1, -2,
         0,  0,  0,  0,  0,  0, -2, -1,  4, -1,
         0,  0,  0,  0,  0,  0, -1, -2, -1,  4
    };

    ASSERT_DOUBLE_ARRAY_EQ(expected_m, eqs.getA()->getData(), 100);
}
