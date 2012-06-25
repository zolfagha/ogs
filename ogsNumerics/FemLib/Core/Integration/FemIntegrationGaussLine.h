
#pragma once

#include "MathLib/Integration/GaussLegendre.h"
#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"

namespace FemLib
{

class FemIntegrationGaussLine : public AbstractFemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        x[0] = MathLib::GaussLegendre::getPoint(getSamplingLevel(), igp);
    }

    double getWeight(size_t igp) const
    {
        return MathLib::GaussLegendre::getWeight(getSamplingLevel(), igp);
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement&, size_t n_sampl_level) const
    {
        return n_sampl_level;
    }
};

}
