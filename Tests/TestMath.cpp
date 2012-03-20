
#include <gtest/gtest.h>

#include <valarray>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"

#include "TestUtil.h"

using namespace MathLib;

template <class T_D0, class T_D1>
class IDifferentiableFunction
{
public:
	  virtual void eval0(T_D0 &x, T_D0 &r) = 0;
	  virtual void eval1(T_D0 &x, T_D1 &j) = 0;
};

class NR1 : public IDifferentiableFunction<double,double>
{
public:
	  void eval0(double &x, double &r) { r = x*x-4.; }
	  void eval1(double &x, double &j) { j = 2*x; }
};

class NR2 : public IDifferentiableFunction<std::valarray<double>,MathLib::Matrix<double> >
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


template<class T_VALUE, class T_JACOBIAN>
class MyNR
{
public:
    template<class F_PROBLEM>
    int solve(F_PROBLEM &fun, T_VALUE &x0, T_VALUE &x_new, size_t max_itr_count=100, double err=1.e-3)
    {
    	T_JACOBIAN jacobian;
    	create_jac(x0, jacobian);
    	T_VALUE r = x0;
    	T_VALUE dx = x0;
    	x_new = x0;

    	bool converged = false;
    	std::cout << "Nonlinear iteration started!" << std::endl;
    	for (size_t i=0; i<max_itr_count; i++) {
        	fun.eval0(x_new, r);
        	fun.eval1(x_new, jacobian);
            if (!check_jac(jacobian)) {
            	std::cout << "->***Warning: Jacobian was evaluated as zero." << std::endl;
            	return -1;
            }
            calculate_dx(r, jacobian, dx);
            x_new += dx;
            printout(i, x_new, dx);
            if (fabs(get_error(dx)) < err) {
                converged = true;
                break;
            }
        }

    	if (converged) return 0;
    	std::cout << "->*** Warning: the iterations didn't converge." << std::endl;

        return -1; //not converged
    }

private:
    inline void create_jac(T_VALUE &x0, T_JACOBIAN &j) {};
    inline bool check_jac(T_JACOBIAN &j);
    inline void calculate_dx(T_VALUE &r, T_JACOBIAN &jac, T_VALUE &dx);
    inline void printout(size_t i, T_VALUE& x_new, T_VALUE& dx);
    inline double get_error(T_VALUE &dx);
};

template<>
inline double MyNR<double, double>::get_error(double& dx) {	return dx; }

template<>
inline void MyNR<double, double>::printout(size_t i, double& x_new, double& dx)
{
	std::cout << "-> " << i <<": x=" << x_new << ", dx=" << dx << std::endl;
}

template<>
inline bool MyNR<double, double>::check_jac(double &j)
{
	return j!=.0;
}

template<>
inline void MyNR<double, double>::calculate_dx(double &r, double &jac, double &dx)
{
    dx = -r/jac;
}

template<>
inline void MyNR<std::valarray<double>, MathLib::Matrix<double> >::create_jac(std::valarray<double> &x0, MathLib::Matrix<double> &j)
{
	j.resize(x0.size(), x0.size());
}

template<>
inline double MyNR<std::valarray<double>, MathLib::Matrix<double> >::get_error(std::valarray<double>& dx)
{
	return dx.max();
}

template<>
inline void MyNR<std::valarray<double>, MathLib::Matrix<double> >::printout(size_t i, std::valarray<double>& x_new, std::valarray<double>& dx)
{
	std::cout << "-> " << i <<": x=(";
	for (size_t i=0; i<x_new.size(); i++) std::cout << x_new[i] << " ";
	std::cout << "), dx=(";
	for (size_t i=0; i<dx.size(); i++) std::cout << dx[i] << " ";
	std::cout << ")" << std::endl;
}

template<>
inline bool MyNR<std::valarray<double>, MathLib::Matrix<double> >::check_jac(MathLib::Matrix<double> &jac)
{
	for (size_t i=0; i<jac.getNRows(); i++)
		for (size_t j=0; j<jac.getNCols(); j++)
			if (jac(i,j)!=.0) return true;
	return false;
}

template<>
inline void MyNR<std::valarray<double>, MathLib::Matrix<double> >::calculate_dx(std::valarray<double> &r, MathLib::Matrix<double> &jac, std::valarray<double> &dx)
{
	MathLib::DenseLinearEquations solver;
	solver.create(r.size());
	for (size_t i=0; i<jac.getNRows(); i++)
		for (size_t j=0; j<jac.getNCols(); j++)
			solver.setA(i,j,jac(i,j));
	for (size_t i=0; i<r.size(); ++i)
		solver.setRHS(i, -1.*r[i]);
	solver.solve();
	double *x = solver.getX();
	for (size_t i=0; i<dx.size(); ++i)
		dx[i] = x[i];
}

TEST(Math, NonlinearNR1)
{
	NR1 f1;
	double x0 = 6.0;
	double x = .0;
	MyNR<double,double> nr;
	nr.solve(f1, x0, x);
	ASSERT_NEAR(2.0, x, 1e-5);
}

TEST(Math, NonlinearNR2)
{
	NR2 f2;
	std::valarray<double> x0(6.0, 2);
	std::valarray<double> x(0.0, 2);
	MyNR<std::valarray<double>,MathLib::Matrix<double> > nr;
	nr.solve(f2, x0, x);

	double my_expect[] = {1., 1.};
	ASSERT_DOUBLE_ARRAY_EQ(my_expect, x, 2, 1e-5);
}

template<class T_VECTOR>
class FResid
{
public:
	void eval(T_VECTOR &x, T_VECTOR &r)
	{
	}
};

template<class T_VECTOR, class T_MATRIX>
class FJacob
{
public:
	void eval(T_VECTOR &x, T_MATRIX &j)
	{
	}
};

template<class T_VECTOR, class T_MATRIX, class T_LINEAR_EQUATION>
class LinearSolverForNR
{
public:
	LinearSolverForNR(T_LINEAR_EQUATION &linear_eqs)
	{
		_linear_equation = &linear_eqs;
	}

    void solve(T_MATRIX &jacobian, T_VECTOR &r, T_VECTOR &dx)
    {
    	_linear_equation->solve();
    }

    T_MATRIX* getJacobian()
    {
    	return _linear_equation->getA();
    }

    T_VECTOR* getResidual()
    {
    	return _linear_equation->getRHS();
    }

private:
    T_LINEAR_EQUATION* _linear_equation;
};


//TEST(Math, NonlinearNR1)
//{
//	typedef std::valarray<double> MyVector;
//	typedef Matrix<double> MyMatrix;
//	typedef FResid<MyVector> MyResid;
//	typedef FJacob<MyVector, MyMatrix> MyJacob;
//	typedef LinearSolverForNR<MyVector, MyMatrix, DenseLinearEquations> MyLinearSolver;
//	MyResid f_r;
//	MyJacob f_j;
//	MyVector x0;
//	MyVector x_new;
//
//
//	NewtoRaphson nr;
//	nr.solve<MyResid, MyJacob, MyMatrix, MyVector, MyLinearSolver, ConvergenceChecker>(f_r, f_j, x0, x_new);
//}

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
