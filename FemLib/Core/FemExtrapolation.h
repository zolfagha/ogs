
#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

#include "MathLib/Vector.h"
#include "MathLib/Integration/GaussLegendre.h"

#include "MeshLib/Core/IElement.h"
#include "FemLib/Core/Integration.h"

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

/**
 * \brief Extrapolation of Gauss point values to nodal values
 */
class FeExtrapolationGaussLinear : public IFeExtrapolationMethod
{
public:
    void extrapolate(IFiniteElement &fe, const std::vector<MathLib::Vector> &gp_values, std::vector<MathLib::Vector> &nodal_values);

    template<typename Tvalue>
    void extrapolate(IFiniteElement &fe, const std::vector<Tvalue> &gp_values, std::vector<Tvalue> &nodal_values);

private:
    double calcXi_p(const MeshLib::IElement& e, FemIntegrationGaussBase& gauss);

    int getLocalIndex(const MeshLib::IElement& e, FemIntegrationGaussBase& gauss, size_t igp);

    void getExtropoGaussPoints(const MeshLib::IElement &e, const int i, double Xi_p, double* unit);
};

} // end namespace
