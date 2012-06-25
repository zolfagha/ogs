
#pragma once

#include "Convergence.h"

namespace MathLib
{

/**
 * \brief Picard method (Fixed-pint iteration method)
 */
class PicardMethod
{
public:
	/// general solver
	/// @tparam F_PROBLEM		Fixed point function class
	/// @tparam T_VALUE			Variable data class
	/// @tparam T_CONVERGENCE	Convergence checker class
    /// @param fun				Fixed point function
    /// @param x0				initial guess
    /// @param x_new			solution
    /// @param x_old			temporal data
    /// @param dx				temporal data
    /// @param error			error tolerance
    /// @param max_itr_count	maximum iteration counts
    template<class F_PROBLEM, class T_VALUE, class T_CONVERGENCE>
    int solve(F_PROBLEM &fun,  const T_VALUE &x0, T_VALUE &x_new, T_VALUE &x_old, T_VALUE &dx, size_t max_itr_count=100, T_CONVERGENCE* convergence=0)
    {
        T_CONVERGENCE _default_convergence;
    	if (convergence==0) convergence = &_default_convergence;
    	x_old = x0;

    	bool converged = false;
    	std::cout << "Picard iteration started!" << std::endl;
    	for (size_t i=0; i<max_itr_count; i++) {
        	fun.eval(x_old, x_new);
        	dx = x_new;
        	dx -= x_old;
            printout(i, x_new, dx);
            if (convergence->check(0, &dx, &x_new)) {
                converged = true;
                break;
            }
            x_old = x_new;
        }

    	if (converged) return 0;
    	std::cout << "->*** Warning: the iterations didn't converge." << std::endl;

        return -1; //not converged
    }

    /// solve scalar problems
    template<class F_PROBLEM>
    int solve(F_PROBLEM &fun, const double x0, double &x_new, double error=1e-6, size_t max_itr_count=100)
    {
    	double x_old, dx;
    	NRCheckConvergence<double,NRErrorNorm1DX> check(error);
    	return solve(fun, x0, x_new, x_old, dx, max_itr_count, &check);
    }

    /// solve vector problems using a direct linear solver
    template<class F_PROBLEM, class T_V, class T_CONVERGENCE>
    int solve(F_PROBLEM &fun, const T_V &x0, T_V &x_new, T_CONVERGENCE* check_error=0, size_t max_itr_count=100)
    {
    	const size_t n = x0.size();
    	T_V x_old(n), dx(n);
    	return solve(fun, x0, x_new, x_old, dx, max_itr_count, check_error);
    }

    /// solve vector problems using a direct linear solver
    template<class F_PROBLEM, class T_V>
    int solve(F_PROBLEM &fun, const T_V &x0, T_V &x_new, double error=1e-6, size_t max_itr_count=100)
    {
    	NRCheckConvergence<T_V,NRErrorNorm1DX> check(error);
    	return solve(fun, x0, x_new, &check, max_itr_count);
    }


private:

    /// print out for debugging
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
inline void PicardMethod::printout(size_t i, double& x_new, double& dx)
{
	std::cout << "-> " << i <<": x=" << x_new << ", dx=" << dx << std::endl;
}

} //end
