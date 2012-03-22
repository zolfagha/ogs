
#include <gtest/gtest.h>

#include <valarray>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"
#include "MathLib/Nonlinear/Picard.h"
#include "MathLib/Nonlinear/BisectionMethod.h"
#include "MathLib/Function/IFunction.h"
#include "TestUtil.h"

using namespace MathLib;

template<class T_F, class T_DF>
class TemplateFixedPointFunctonUsingC1Equation
{
	MathLib::DenseLinearEquations solver;
	std::valarray<double> b;
	std::valarray<double> tmp;
	T_F nr;
	T_DF nr_df;
public:
	TemplateFixedPointFunctonUsingC1Equation(size_t n)
	{
		solver.create(n);
		b.resize(n);
		tmp.resize(n);
	}
	// x = g(x) = x - f'(x)^-1 * f(x)
	// f'*g = f'*x - f
	void eval(std::valarray<double> &x, std::valarray<double> &x_new)
	{
		MathLib::Matrix<double> *j = solver.getA();
		solver.reset();
		b = 0;
		tmp = 0;

		nr_df.eval(x, *j);
		nr.eval(x, b);

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



// Example problems (one variable) ------------------------------------------------------------------
// f(x) = x*x -4 = 0
// x = 2,-2
class NL1_NR : public AbstractDefaultCloneFunction<double, double, NL1_NR>
{
public:
	// r = f(x)
	void eval(const double &x, double &r) { r = x*x-4.; }
};

class NL1_NR_D1 : public AbstractDefaultCloneFunction<double, double, NL1_NR_D1>
{
public:
	// j = f'(x)
	void eval(const double &x, double &j) { j = 2*x; }
};

class NL1_PC : public AbstractDefaultCloneFunction<double, double, NL1_PC>
{
public:
	// x = g(x) = x - 1/f'(x) * f(x)
	//void eval0(double &x, double &x_new) { x_new = 0.5*(2./x+x); }
	void eval(const double &x, double &x_new) { x_new = (x*x+4.)/(2.*x); }
};


// Example problems (2 variables) ------------------------------------------------------------------
// 3x-y=-2
// 2x^2-y=0
// (x,y) = (-1/2, 1/2) and (2, 8)
class NL2_NR
{
public:
	void eval(std::valarray<double> &x, std::valarray<double> &r)
	{
		r[0] = 3*x[0]-x[1]+2.;
		r[1] = 2*x[0]*x[0]-x[1];
	}
};

class NL2_NR_D1
{
public:
	void eval(std::valarray<double> &x, MathLib::Matrix<double> &j)
	{
		j(0,0) = 3.;
		j(0,1) = -1.0;
		j(1,0) = 4.*x[0];
		j(1,1) = -1.0;
	}
};

typedef TemplateFixedPointFunctonUsingC1Equation<NL2_NR,NL2_NR_D1> NL2_PC1;


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
	void eval(std::valarray<double> &x, std::valarray<double> &x_new)
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


class NL3_NR
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
	  _linear_solver.getOption().ls_method = MathLib::LIS_option::BiCGSTAB;
  }

  void eval(std::valarray<double> &x, std::valarray<double> &r)
  {
		r[0] = 3*x[0]-x[1]+2.;
		r[1] = 2*x[0]*x[0]-x[1];
  }

  MathLib::CRSLisSolver* getLinearSolver() {return &_linear_solver;};


private:
  MathLib::CRSLisSolver _linear_solver;
};

class NL3_NR_D1
{
public:
	NL3_NR_D1(MathLib::CRSLisSolver* solver) : _linear_solver(solver)
	{
	}

	void eval(std::valarray<double> &x, MathLib::CRSMatrix<double, signed> &j)
	{
		j.setValue(0, 0, 3.);
		j.setValue(0, 1, -1.0);
		j.setValue(1, 0, 4.*x[0]);
		j.setValue(1, 1,  -1.0);
	}

	MathLib::CRSLisSolver* getLinearSolver() {return _linear_solver;};


private:
  MathLib::CRSLisSolver* _linear_solver;
};
//------------------------------------------------------------------

TEST(Math, NonlinearNR_double)
{
	NL1_NR f1;
	NL1_NR_D1 df1;
	double x0 = 6.0;
	double x = .0;
	NewtonRaphsonMethod nr;
	nr.solve(f1, df1, x0, x);
	ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearNR_dense)
{
	typedef std::valarray<double> MyVector;

	NL2_NR f2;
	NL2_NR_D1 df2;
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	NewtonRaphsonMethod nr;
	nr.solve(f2,df2, x0, x);

	double my_expect[] = {2., 8.};
//	double my_expect[] = {1., 1.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearNR_sparse)
{
	typedef std::valarray<double> MyVector;
	typedef MathLib::CRSMatrix<double, signed> MyMatrix;
	typedef NRCheckConvergence<MyVector,NRErrorNorm1DX > MyConverge;
	typedef NewtonDxSolverVector<MyVector,MyMatrix> MyDxSolver;

	NL3_NR f;
	NL3_NR_D1 df(f.getLinearSolver());
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	MyVector r(2), dx(2);
	MyDxSolver dx_solver(f.getLinearSolver());
	MyMatrix* j = f.getLinearSolver()->getA();

	NewtonRaphsonMethod nr;
	nr.solve<NL3_NR,NL3_NR_D1,MyVector,MyMatrix,MyDxSolver,MyConverge>(f, &df, x0, x, r, *j, dx, dx_solver);

	double my_expect[] = {2., 8.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearPicard_double)
{
	NL1_PC f1;
	double x0 = 6.0;
	double x = .0;
	PicardMethod nr;
	nr.solve(f1, x0, x);
//	ASSERT_NEAR(sqrt(2.0), x, 1e-5);
	ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearPicard_dense1)
{
	typedef std::valarray<double> MyVector;

	NL2_PC1 f2(2);
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	PicardMethod nr;
	nr.solve(f2, x0, x);

//	double my_expect[] = {1., 1.};
	double my_expect[] = {2., 8.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearPicard_dense2)
{
	typedef std::valarray<double> MyVector;

	NL2_PC2 f2(2);
	MyVector x0(6.0, 2);
	MyVector x(0.0, 2);
	PicardMethod nr;
	nr.solve(f2, x0, x);

	double my_expect[] = {-0.5, 0.5};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
//    ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

TEST(Math, NonlinearBS_double)
{
	NL1_NR f1;
	double x = .0;
	double a = .0;
	double c = 6.0;
	BisectionMethod nr;
	nr.solve(f1, a, c, x);
	ASSERT_NEAR(2.0, x, 1e-5);
}

