
#pragma once

#include "FemFunction.h"

namespace FemLib
{
template<typename Tvalue, typename Tpos>
class IFemExtrapolation
{
public:
    virtual void extrapolate(FEMIntegrationPointFunction<Tvalue, Tpos> &ele, FEMNodalFunction<Tvalue, Tpos> &nod) = 0;
};

template<typename Tvalue, typename Tpos>
class FEMExtrapolationAverage : public IFemExtrapolation<typename Tvalue, typename Tpos>
{
public:
    void extrapolate(FEMIntegrationPointFunction<Tvalue, Tpos> &ele, FEMNodalFunction<Tvalue, Tpos> &nod)
    {
        throw std::exception("The method or operation is not implemented.");
    }
};

template<typename Tvalue, typename Tpos>
class FEMExtrapolationLinear : public IFemExtrapolation<typename Tvalue, typename Tpos>
{
    void extrapolate(FEMIntegrationPointFunction<Tvalue, Tpos> &ele, FEMNodalFunction<Tvalue, Tpos> &nod)
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

template<typename Tvalue, typename Tpos>
class FEMExtrapolationFactory
{
public:
    static IFemExtrapolation<Tvalue, Tpos>* create(FEMExtrapolationMethod::type tp)
    {
        switch (tp) {
            case FEMExtrapolationMethod::Linear:
                return new FEMExtrapolationLinear<Tvalue, Tpos>();
            case FEMExtrapolationMethod::Average:
                return new FEMExtrapolationAverage<Tvalue, Tpos>();
        }
        return 0;
    };
};


}
