
#pragma once

#include <vector>

#include "MathLib/Vector.h"

namespace FemLib
{

class IFiniteElement;

/**
 * \brief Extrapolation of integration point values to nodal values
 */
class IFeExtrapolationMethod
{
public:
	virtual ~IFeExtrapolationMethod() {};
    virtual void extrapolate(IFiniteElement &fe, const std::vector<MathLib::Vector> &gp_values, std::vector<MathLib::Vector> &nodal_values) = 0;
};

} // end
