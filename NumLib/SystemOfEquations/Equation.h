
#pragma once

#include <vector>
#include <valarray>

#include "Base/CodingTools.h"
#include "MathLib/LinAlg/Dense/Matrix.h"

#include "Variable.h"


namespace NumLib
{

struct TimeODEComponent
{
	typedef Matrix<double> MatrixType;
	typedef std::valarray<double> VectorType;

	MatrixType* mass;
	MatrixType* stiffness;
	VectorType* source;
};

struct TimeODE
{
	typedef Matrix<double> MatrixType;
	typedef std::valarray<double> VectorType;

	std::vector<Variable> variables;
	std::vector<TimeODEComponent> componets;
	Variable primary_variable;
	size_t n_dof;

	void addVariable(Variable &var, TimeODEComponent &comp, bool is_primary=false)
	{
		variables.push_back(var);
		componets.push_back(comp);
		if (is_primary) {
			primary_variable = var;
			n_dof = var.n_dof;
		}
	}

};

struct LinearComponent
{
	typedef Matrix<double> MatrixType;
	typedef std::valarray<double> VectorType;

	MatrixType* stiffness;
	VectorType* source;

	LinearComponent(MatrixType* m, VectorType* v = 0) : stiffness(m), source(v) {};
};

struct LinearEquation
{
	std::vector<Variable*> variables;
	std::vector<LinearComponent*> componets;
	Variable* primary_variable;
	size_t n_dof;

	void addVariable(Variable &var, LinearComponent &comp, bool is_primary=false)
	{
		variables.push_back(&var);
		componets.push_back(&comp);
		if (is_primary) {
			primary_variable = &var;
			n_dof = var.n_dof;
		}
	}

    virtual ~LinearEquation()
    {
        Base::releaseObjectsInStdVector(componets);
    }
};

} //end

