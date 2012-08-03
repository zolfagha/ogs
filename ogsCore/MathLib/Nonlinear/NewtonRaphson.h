/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NewtonRaphson.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include "logog/include/logog.hpp"

#include "MathLib/LinAlg/MatrixBase.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "Convergence.h"

namespace MathLib
{


template<class F_JACOBIAN, class T_LINEAR_SOLVER>
class NewtonFunctionDXVector
{
private:
    F_JACOBIAN* _f_j;
    T_LINEAR_SOLVER* _linear_solver;
    typedef typename T_LINEAR_SOLVER::MatrixType JacobianType;
public:
    NewtonFunctionDXVector(F_JACOBIAN &f, T_LINEAR_SOLVER &linear_solver)
        : _f_j(&f), _linear_solver(&linear_solver) {};

    template<class T_VALUE>
    void eval(const T_VALUE &x, const T_VALUE &r, T_VALUE &dx)
    {
        const size_t n_rows = _linear_solver->getDimension();
        JacobianType* jac = _linear_solver->getA();
        _f_j->eval(x, *jac);
        //if (!check_jac(jac)) return false;

        // rhs = -r
        for (size_t i=0; i<n_rows; ++i)
            _linear_solver->setRHS(i, -1.*r[i]);

        // solve J dx = -r
        _linear_solver->solve();

        // get dx
        double *u = _linear_solver->getX();
        for (size_t i=0; i<n_rows; ++i)
            dx[i] = u[i];
    }

private:
    inline bool check_jac(JacobianType &jac)
    {
        const size_t n = _linear_solver->getDimension();
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                if (jac(i,j)!=.0) return true;
        return false;
    }
};

template<class F_JACOBIAN>
class NewtonFunctionDXScalar
{
private:
    F_JACOBIAN* _f_j;
public:
    NewtonFunctionDXScalar(F_JACOBIAN &f) : _f_j(&f) {};

    void eval(const double &x, const double &r, double &dx)
    {
        double j;
        _f_j->eval(x, j);
        if (j==0) return; //TODO error
        dx = -r / j;
    }
};


/**
 * \brief Newton-Raphson method
 */
class NewtonRaphsonMethod
{
public:
    /// general solver
    /// @tparam F_RESIDUALS
    /// @tparam F_DX
    /// @tparam T_VALUE
    /// @tparam T_CONVERGENCE
    template<class F_RESIDUALS, class F_DX, class T_VALUE, class T_CONVERGENCE>
    int solve(F_RESIDUALS &f_residuals, F_DX &f_dx, const T_VALUE &x0, T_VALUE &x_new, T_VALUE &r, T_VALUE &dx, size_t max_itr_count=100, T_CONVERGENCE* convergence=0)
    {
        T_CONVERGENCE _default_convergence;
        if (convergence==0) convergence = &_default_convergence;
        r = .0;
        x_new = x0;

        INFO("Newton-Raphson iteration started!");
        bool converged = false;
        size_t itr_cnt=0;
        f_residuals.eval(x_new, r);
        //printout(0, x_new, r, dx);
        if (!convergence->check(&r, &dx, &x_new)) {
            for (itr_cnt=0; itr_cnt<max_itr_count; itr_cnt++) {
                f_dx.eval(x_new, r, dx);
                x_new += dx;
                f_residuals.eval(x_new, r);
                //printout(itr_cnt, x_new, r, dx);
                if (convergence->check(&r, &dx, &x_new)) {
                    converged = true;
                    break;
                }
            }
        }

        INFO("------------------------------------------------------------------");
        INFO("*** Newton-Raphson nonlinear solver computation result");
        if (max_itr_count==1) {
            INFO("status    : iteration not required");
        } else {
            INFO("status    : %s", (converged ? "converged" : "***ERROR - DID NOT CONVERGE!"));
        }
        INFO("iteration : %d/%d", itr_cnt, max_itr_count);
        INFO("residuals : %1.3e (tolerance=%1.3e)", convergence->getError(), convergence->getTolerance());
        INFO("------------------------------------------------------------------");

        if (converged) return 0;
        std::cout << "->*** Warning: the iterations didn't converge." << std::endl;

        return -1; //not converged
    }

    /// solve scalar problems
    template<class F_RESIDUALS, class F_JACOBIAN>
    int solve(F_RESIDUALS &f_residuals, F_JACOBIAN &f_jac, const double &x0, double &x_new, double error=1e-6, size_t max_itr_count=100)
    {
        NewtonFunctionDXScalar<F_JACOBIAN> f_dx(f_jac);
        double r, dx;
        NRCheckConvergence<double,NRErrorAbsResMNormOrRelDxMNorm> check(error);
        return solve(f_residuals, f_dx, x0, x_new, r, dx, max_itr_count, &check);
    }

    /// solve vector problems using a direct linear solver
    template<class F_RESIDUALS, class F_JACOBIAN, class T_V, class T_CONVERGENCE>
    int solve(F_RESIDUALS &f_residuals, F_JACOBIAN &f_jac, const T_V &x0, T_V &x_new, T_CONVERGENCE* check_error=0, size_t max_itr_count=100)
    {
        const size_t n = x0.size();
        T_V r(n), dx(n);
        MathLib::DenseLinearEquations dense;
        dense.create(n);
        NewtonFunctionDXVector<F_JACOBIAN, MathLib::DenseLinearEquations> f_dx(f_jac, dense);
        return solve(f_residuals, f_dx, x0, x_new, r, dx, max_itr_count, check_error);
    }

    /// solve vector problems using a direct linear solver
    template<class F_RESIDUALS, class F_JACOBIAN, class T_V>
    int solve(F_RESIDUALS &f_residuals, F_JACOBIAN &f_jac, const T_V &x0, T_V &x_new, double error=1e-6, size_t max_itr_count=100)
    {
        NRCheckConvergence<T_V,NRErrorAbsResMNormOrRelDxMNorm> check(error);
        return solve(f_residuals, f_jac, x0, x_new, &check, max_itr_count);
    }

private:
    template<class T_VALUE>
    inline void printout(size_t i, T_VALUE &x_new, T_VALUE &r, T_VALUE &dx)
    {
        std::cout << "-> " << i <<": ";
        std::cout << "r=(";
        for (size_t i=0; i<dx.size(); i++) std::cout << r[i] << " ";
#if 0
        std::cout << "-> " << i <<": x=(";
        for (size_t i=0; i<x_new.size(); i++) std::cout << x_new[i] << " ";
        std::cout << "), r=(";
        for (size_t i=0; i<dx.size(); i++) std::cout << r[i] << " ";
        std::cout << "), dx=(";
        for (size_t i=0; i<dx.size(); i++) std::cout << dx[i] << " ";
        std::cout << ")" << std::endl;
#endif
    }
};

#if 0
// template specialization
template<>
inline void NewtonRaphsonMethod::printout(size_t i, double &x_new, double &r, double &dx)
{
    std::cout << "-> " << i << ": x=" << x_new << ", r=" << r << ", dx=" << dx << std::endl;
}
#endif

} //end

