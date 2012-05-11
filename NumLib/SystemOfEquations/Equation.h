
#pragma once

#include <vector>
#include <valarray>

#include "Base/CodingTools.h"
#include "NumLib/DataType.h"
#include "Variable.h"


namespace NumLib
{

struct TimeODEComponent
{
	LocalMatrix* mass;
	LocalMatrix* stiffness;
	LocalVector* source;
};

struct TimeODE
{
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
	LocalMatrix* stiffness;
	LocalVector* source;

	LinearComponent(LocalMatrix* m, LocalVector* v = 0) : stiffness(m), source(v) {};
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

