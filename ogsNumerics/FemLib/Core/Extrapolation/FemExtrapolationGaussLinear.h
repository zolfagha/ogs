
#pragma once

#include <vector>

#include "MeshLib/Core/IElement.h"

#include "FemLib/Core/DataType.h"
#include "FemLib/Core/Integration/Integration.h"
#include "IFeExtrapolationMethod.h"

namespace FemLib
{

/**
 * \brief Extrapolation of Gauss point values to nodal values
 */
class FeExtrapolationGaussLinear : public IFeExtrapolationMethod
{
public:
    void extrapolate(IFiniteElement &fe, const std::vector<LocalVector> &gp_values, std::vector<LocalVector> &nodal_values);

    template<typename Tvalue>
    void extrapolate(IFiniteElement &fe, const std::vector<Tvalue> &gp_values, std::vector<Tvalue> &nodal_values);

private:
    double calcXi_p(const MeshLib::IElement& e, AbstractFemIntegrationGaussBase& gauss);

    int getLocalIndex(const MeshLib::IElement& e, AbstractFemIntegrationGaussBase& gauss, size_t igp);

    void getExtropoGaussPoints(const MeshLib::IElement &e, const int i, double Xi_p, double* unit);
};

} // end namespace
