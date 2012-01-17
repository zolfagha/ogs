
#pragma once

#include <vector>
#include "MeshLib/Core/IMesh.h"

namespace FemLib
{
class IFemIntegration
{
public:
//    void integrate(int integrant, double val);
};

class FemIntegrationGauss : IFemIntegration
{
public:
    FemIntegrationGauss() {};
    FemIntegrationGauss(MeshLib::ElementType ele_type, int sampling_level) {
    }

    size_t getNumberOfSamplingPoints();
    std::vector<int> getSamplingPoints();
    double* getSamplingPoint(int igp);
    double getWeight(int igp);
};

}
