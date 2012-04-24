
#pragma once

#include "DiscreteLib/Core/DiscreteVector.h"
#include "NonlinearSolverOption.h"

namespace NumLib
{

/**
 * \brief Interface to nonlinear solvers
 */
class INonlinearSolver
{
public:
	typedef DiscreteLib::DiscreteVector<double> VectorType;

	virtual ~INonlinearSolver() {};

	virtual void solve(const VectorType &x_0, VectorType &x_new) = 0;

	NonlinerSolverOption& getOption() {return _option;};
	void setOption(const NonlinerSolverOption &option) {_option = option;};
protected:
	NonlinerSolverOption _option;
};

} //end

