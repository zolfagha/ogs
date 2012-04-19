
#pragma once

#include <vector>

#include "MathLib/Vector.h"

#include "MeshLib/Core/IElement.h"

#include "IFeExtrapolationMethod.h"

namespace FemLib
{
/**
 * \brief Extrapolation of integration point values by taking averages
 */
class FeExtrapolationAverage : public IFeExtrapolationMethod
{
public:
    void extrapolate(IFiniteElement &fe, const std::vector<MathLib::Vector> &gp_values, std::vector<MathLib::Vector> &nodal_values);
    template<typename Tvalue>
    void extrapolate(IFiniteElement &fe, const std::vector<Tvalue> &gp_values, std::vector<Tvalue> &nodal_values);
};

} // end namespace
