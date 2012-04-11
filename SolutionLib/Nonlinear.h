
#pragma once

#include "MathLib/Nonlinear/Picard.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"

namespace SolutionLib
{

/**
 * \brief Linear
 */
template <class T_LINEAR_FUNCTION>
class Linear
{
	T_LINEAR_FUNCTION* _linear_f;
public:
	Linear(T_LINEAR_FUNCTION &linear_f) : _linear_f(&linear_f) {};

	template <class T_VALUE>
	void solve(T_VALUE &x_0, T_VALUE &x_new)
	{
		_linear_f->eval(x_0, x_new);
	}
};

/**
 * \brief Picard
 */
template <class T_LINEAR_FUNCTION>
class Picard
{
	T_LINEAR_FUNCTION* _linear_f;
public:
	Picard(T_LINEAR_FUNCTION &linear_f) : _linear_f(&linear_f) {};

	template <class T_VALUE>
	void solve(T_VALUE &x_0, T_VALUE &x_new)
	{
        MathLib::PicardMethod picard;
        picard.solve(*_linear_f, x_0, x_new);
	}
};


template <class T_LINEAR_FUNCTION>
class NewtonRaphson
{
	T_LINEAR_FUNCTION* _f_J_r;
public:
	NewtonRaphson(T_LINEAR_FUNCTION &f_J_r) : _f_J_r(&f_J_r) {};

	template <class T_VALUE>
	void solve(T_VALUE &x_0, T_VALUE &x_new)
	{
        //MathLib::NewtonRaphsonMethod nr;
        //nr.solve(_f_r, _f_dx, x_0, x_new);
	}
};

//template <class F_RESIDUALS, class F_DX>
//class NewtonRaphson
//{
//	F_RESIDUALS* _f_r;
//	F_DX* _f_dx;
//public:
//	NewtonRaphson(F_RESIDUALS &f_r, F_DX &f_dx) : _f_r(&f_r), _f_dx(f_dx) {};
//
//	template <class T_VALUE>
//	void solve(T_VALUE &x_0, T_VALUE &x_new)
//	{
//        MathLib::NewtonRaphsonMethod nr;
//        nr.solve(_f_r, _f_dx, x_0, x_new);
//	}
//};

} //end

