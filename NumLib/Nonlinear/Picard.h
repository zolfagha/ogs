
#pragma once

#include "MathLib/Nonlinear/Picard.h"
#include "DiscreteLib/Core/DiscreteVector.h"
#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Picard
 */
template <class F_LINEAR>
class Picard : public INonlinearSolver
{
	F_LINEAR* _linear_f;
public:
	explicit Picard(F_LINEAR* linear_f) : _linear_f(linear_f) {};
	virtual ~Picard() {};

	virtual void solve(const VectorType &x_0, VectorType &x_new)
	{
        MathLib::PicardMethod picard;
        picard.solve(*_linear_f, x_0, x_new);
	}
};

} //end

