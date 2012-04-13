
#pragma once

#include <vector>

#include "Base/CodingTools.h"

#include "Variable.h"
#include "Equation.h"


namespace MathLib
{

class SystemOfEquations
{
public:
	SystemOfEquations() {};

	virtual ~SystemOfEquations() 
    {
        Base::releaseObjectsInStdVector(_equations);
        Base::releaseObjectsInStdVector(_variables);
    };

	void addEquation(LinearEquation &eq)
	{
		_equations.push_back(&eq);
		for (size_t i=0; i<eq.variables.size(); i++) {
			if (!this->hasVariable(*eq.variables[i])) {
				this->_variables.push_back(eq.variables[i]);
			}
		}
	}

	size_t getNumberOfEquations() const { return _equations.size(); };

	LinearEquation* getEquation(size_t i) {return _equations[i];};

	std::vector<LinearEquation*>* getListOfEquations() {return &_equations;};

	size_t getNumberOfVariables() const { return _variables.size(); };

	Variable* getVariable(size_t i) {return _variables[i];};

	std::vector<Variable*>* getListOfVariables() {return &_variables;};

	bool hasVariable(const Variable &var) const
	{
		for (size_t i=0; i<_variables.size(); i++) {
			if (var.id == _variables[i]->id) return true;
		}
		return false;
	}

	bool isReadyToSolve() const
	{
		return _variables.size()==_equations.size();
	};

private:
	std::vector<Variable*> _variables;
	std::vector<LinearEquation*> _equations;
};

} //end
