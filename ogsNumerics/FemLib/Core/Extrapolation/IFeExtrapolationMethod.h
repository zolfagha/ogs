
#pragma once

#include <vector>

#include "FemLib/Core/DataType.h"

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
    virtual void extrapolate(IFiniteElement &fe, const std::vector<LocalVector> &gp_values, std::vector<LocalVector> &nodal_values) = 0;
};

} // end
