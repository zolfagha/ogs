
#pragma once

#include <string>

#include "BaseLib/Options.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteVector.h"
#include "NonlinearSolver.h"
#include "NonlinearSolverOption.h"

namespace NumLib
{

template <class F_LINEAR, class F_R, class F_DX>
class TemplateDiscreteNonlinearSolver
{
public:
	typedef DiscreteLib::IDiscreteVector<double> VectorType;

	TemplateDiscreteNonlinearSolver(DiscreteLib::DiscreteSystem* dis_sys, F_LINEAR* f_l, F_R* f_r, F_DX* f_dx)
	: _dis_sys(dis_sys), _f_l(f_l), _f_r(f_r), _f_dx(f_dx), _solver(0)
	{
	}

	virtual ~TemplateDiscreteNonlinearSolver()
	{
		BaseLib::releaseObject(_solver);
	}

    void setOption(const BaseLib::Options &option)
    {
        const BaseLib::Options *op = option.getSubGroup("Nonlinear");
        if (op==0) return;

        if (op->hasOption("solver_type"))
            _option.solver_type = _option.getSolverType(op->getOption("solver_type"));
        if (op->hasOption("error_tolerance"))
            _option.error_tolerance = op->getOption<double>("error_tolerance");
        if (op->hasOption("max_iteration_step"))
            _option.max_iteration = op->getOption<int>("max_iteration_step");
    }

    void setOption(const NonlinerSolverOption &option) { _option = option; }

    NonlinerSolverOption &getOption() const { return _option; }

	void solve(const VectorType &x_0, VectorType &x_new)
	{
		if (_solver==0)
			_solver = create(_option.solver_type);
		_solver->solve(x_0, x_new);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(TemplateDiscreteNonlinearSolver);

	INonlinearSolver* create(NonlinerSolverOption::SolverType type)
	{
		INonlinearSolver* solver = 0;
		switch (type)
		{
		case NonlinerSolverOption::Linear:
			solver = new Linear<F_LINEAR>(_f_l);
			break;
		case NonlinerSolverOption::Picard:
			solver = new Picard<F_LINEAR>(_dis_sys, _f_l);
			break;
		case NonlinerSolverOption::Newton:
			solver = new NewtonRaphson<F_R, F_DX>(_dis_sys, _f_r, _f_dx);
			break;
		default:
			break;
		}
		solver->setOption(_option);
		return solver;
	}


private:
    NonlinerSolverOption _option;
	DiscreteLib::DiscreteSystem* _dis_sys;
    F_LINEAR* _f_l;
    F_R* _f_r;
    F_DX* _f_dx;
    INonlinearSolver* _solver;
};

} //end

