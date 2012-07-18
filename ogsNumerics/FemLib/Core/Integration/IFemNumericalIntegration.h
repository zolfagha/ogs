
#pragma once

#include "MeshLib/Core/IElement.h"

namespace FemLib
{

/**
 * \brief Interface for integration methods used in FEM
 */
class IFemNumericalIntegration
{
public:
    virtual ~IFemNumericalIntegration() {};
    virtual void initialize(MeshLib::IElement &e, size_t n_sampl_level) = 0;
    virtual size_t getNumberOfSamplingPoints() const = 0;
    virtual void getSamplingPoint(size_t igp, double *) const = 0;
    virtual double getWeight(size_t igp) const = 0;
};

}
