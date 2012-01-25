
#pragma once

#include <vector>
#include "MeshLib/Core/IMesh.h"

namespace FemLib
{
class IFemIntegration
{
public:
//    void integrate(int integrant, double val);

    virtual size_t getNumberOfSamplingPoints() const = 0;
    virtual const double* getSamplingPoint(int igp) const = 0;
    virtual double getWeight(int igp) const = 0;
};

class FemIntegrationGauss : public IFemIntegration
{
public:
    FemIntegrationGauss() {};
    FemIntegrationGauss(MeshLib::ElementType ele_type, int sampling_level) {
    }

    size_t getNumberOfSamplingPoints() const;
    std::vector<int> getSamplingPoints() const;
    const double* getSamplingPoint(int igp) const;
    double getWeight(int igp) const;
};

}
