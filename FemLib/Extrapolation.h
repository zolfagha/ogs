
#pragma once

#include "FemFunction.h"

namespace FemLib
{
template<typename Tvalue>
class IFemExtrapolation
{
public:
    virtual void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele, TemplateFEMNodalFunction<Tvalue> &nod) = 0;
};

template<typename Tvalue>
class FEMExtrapolationAverage : public IFemExtrapolation<typename Tvalue>
{
public:
    void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele, TemplateFEMNodalFunction<Tvalue> &nod)
    {
        throw std::exception("The method or operation is not implemented.");
    }
};

template<typename Tvalue>
class FEMExtrapolationLinear : public IFemExtrapolation<typename Tvalue>
{
    void extrapolate(TemplateFEMIntegrationPointFunction<Tvalue> &ele, TemplateFEMNodalFunction<Tvalue> &nod)
    {
        throw std::exception("The method or operation is not implemented.");
    }
};


struct FEMExtrapolationMethod
{
    enum type {
        Linear = 1,
        Average = 2,
        LeastSquare = 3,
        INVALID = -1
    };
};

template<typename Tvalue>
class FEMExtrapolationFactory
{
public:
    static IFemExtrapolation<Tvalue>* create(FEMExtrapolationMethod::type tp)
    {
        switch (tp) {
            case FEMExtrapolationMethod::Linear:
                return new FEMExtrapolationLinear<Tvalue>();
            case FEMExtrapolationMethod::Average:
                return new FEMExtrapolationAverage<Tvalue>();
        }
        return 0;
    };
};


}
