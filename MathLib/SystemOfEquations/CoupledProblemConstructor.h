
#pragma once

#include <vector>

#include "SystemOfEquations.h"
#include "MyCouplingEQS.h"

namespace MathLib
{

class CoupledProblemConstructor
{
public:
	MyCouplingEQS* createPartitionedProblem(SystemOfEquations &sys, const std::vector<Variable*> &active_vars)
	{
		std::vector<LinearEquation*> active_eqs;
		for (size_t i=0; i<active_vars.size(); i++) {
			for (size_t j=0; j<sys.getNumberOfEquations(); j++) {
				if (active_vars[i]->id == sys.getEquation(j)->primary_variable->id)
					active_eqs.push_back(sys.getEquation(j));
			}
		}
		MyCouplingEQS* eqs = new MyCouplingEQS(active_vars, active_eqs);
		return eqs;
	}

	MyCouplingEQS* createMonolithicProblem(SystemOfEquations &sys)
	{
		MyCouplingEQS* eqs = new MyCouplingEQS(*sys.getListOfVariables(), *sys.getListOfEquations());
		return eqs;
	}
};


} //end
