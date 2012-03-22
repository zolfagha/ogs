
#pragma once

#include "MathLib/LinAlg/MatrixBase.h"
#include "Convergence.h"

namespace MathLib
{

/**
 * Solve J dx = -r
 */
class NewtonDxSolverScalar
{
public:
    inline bool solve(double &r, double &jac, double &dx)
    {
    	if (jac==0) return false;
    	dx = -r / jac;
    	return true;
    }
    void reset() {};
};

/**
 * Solve J dx = -r
 */
template <typename T_VEC, typename T_JAC>
class NewtonDxSolverVector
{
private:
	MathLib::ILinearEquations* _linear_solver;

public:
	NewtonDxSolverVector(MathLib::ILinearEquations* linear_solver) : _linear_solver(linear_solver) {};

    void reset()
    {
    	_linear_solver->reset();
    };

    inline bool solve(T_VEC &r,T_JAC &jac, T_VEC &dx)
    {
    	if (!check_jac(jac)) return false;

    	const size_t n_rows = jac.getNRows();
//    	for (size_t i=0; i<n_rows; i++)
//    		for (size_t j=0; j<n_rows; j++)
//    			_linear_solver->setA(i,j,jac(i,j));

    	// rhs = -r
    	for (size_t i=0; i<n_rows; ++i)
    		_linear_solver->setRHS(i, -1.*r[i]);

    	// solve J dx = -r
    	_linear_solver->solve();

    	// get dx
    	double *x = _linear_solver->getX();
    	for (size_t i=0; i<n_rows; ++i)
    		dx[i] = x[i];

    	return true;
    }

    inline bool check_jac(T_JAC &jac)
    {
    	for (size_t i=0; i<jac.getNRows(); i++)
    		for (size_t j=0; j<jac.getNCols(); j++)
    			if (jac(i,j)!=.0) return true;
    	return false;
    }

};

/**
 * \brief Newton-Raphson method
 */
class NewtonRaphsonMethod
{
public:
	/// general solver
    template<class F_PROBLEM, class F_PROBLEM_DIFF, class T_VALUE, class T_JACOBIAN, class T_DX_SOLVER, class T_CONVERGENCE>
    int solve(F_PROBLEM &fun, F_PROBLEM_DIFF* df, T_VALUE &x0, T_VALUE &x_new, T_VALUE &r, T_JACOBIAN &jacobian, T_VALUE &dx, T_DX_SOLVER &dx_solver, size_t max_itr_count=100, T_CONVERGENCE* convergence=0)
    {
        T_CONVERGENCE _default_convergence;
    	if (convergence==0) convergence = &_default_convergence;
    	r = .0;
    	x_new = x0;

    	bool converged = false;
    	std::cout << "Nonlinear iteration started!" << std::endl;
    	for (size_t i=0; i<max_itr_count; i++) {
    		dx_solver.reset();
        	fun.eval(x_new, r);
        	df->eval(x_new, jacobian);
            if (!dx_solver.solve(r, jacobian, dx)) {
            	std::cout << "->***Warning: Jacobian was evaluated as zero." << std::endl;
            	return -1;
            }
            x_new += dx;
            printout(i, x_new, dx);
            if (convergence->check(&r, &dx, &x_new)) {
                converged = true;
                break;
            }
        }

    	if (converged) return 0;
    	std::cout << "->*** Warning: the iterations didn't converge." << std::endl;

        return -1; //not converged
    }

    /// solve scalar problems
    template<class F_PROBLEM, class F_PROBLEM_DIFF>
    int solve(F_PROBLEM &fun, F_PROBLEM_DIFF &df, double &x0, double &x_new, double error=1e-6, size_t max_itr_count=100)
    {
    	double r, j, dx;
    	NewtonDxSolverScalar dx_solver;
    	NRCheckConvergence<double,NRErrorNorm1DX> check(error);
    	return solve(fun, &df, x0, x_new, r, j, dx, dx_solver, max_itr_count, &check);
    }

    /// solve vector problems using a direct linear solver
    template<class F_PROBLEM, class F_PROBLEM_DIFF, class T_V, class T_CONVERGENCE>
    int solve(F_PROBLEM &fun, F_PROBLEM_DIFF &df, T_V &x0, T_V &x_new, T_CONVERGENCE* check_error=0, size_t max_itr_count=100)
    {
    	const size_t n = x0.size();
    	T_V r(n), dx(n);
    	MathLib::DenseLinearEquations dense;
    	dense.create(n);
    	NewtonDxSolverVector<T_V, MathLib::Matrix<double> > dx_solver(&dense);
    	return solve(fun, &df, x0, x_new, r, *dense.getA(), dx, dx_solver, max_itr_count, check_error);
    }

    /// solve vector problems using a direct linear solver
    template<class F_PROBLEM, class F_PROBLEM_DIFF, class T_V>
    int solve(F_PROBLEM &fun, F_PROBLEM_DIFF &df, T_V &x0, T_V &x_new, double error=1e-6, size_t max_itr_count=100)
    {
    	NRCheckConvergence<T_V,NRErrorNorm1DX> check(error);
    	return solve(fun, df, x0, x_new, &check, max_itr_count);
    }

private:
    template<class T_VALUE>
    inline void printout(size_t i, T_VALUE& x_new, T_VALUE& dx)
    {
    	std::cout << "-> " << i <<": x=(";
    	for (size_t i=0; i<x_new.size(); i++) std::cout << x_new[i] << " ";
    	std::cout << "), dx=(";
    	for (size_t i=0; i<dx.size(); i++) std::cout << dx[i] << " ";
    	std::cout << ")" << std::endl;
    }
};

// template specialization
template<>
inline void NewtonRaphsonMethod::printout(size_t i, double& x_new, double& dx)
{
	std::cout << "-> " << i <<": x=" << x_new << ", dx=" << dx << std::endl;
}

} //end

