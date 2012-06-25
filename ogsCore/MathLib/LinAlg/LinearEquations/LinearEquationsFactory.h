
#pragma once

#include "ILinearEquations.h"
#include "DenseLinearEquations.h"
#include "SparseLinearEquations.h"
#include "LisInterface.h"

namespace MathLib
{

struct LinearEquationsType
{
    enum type {
        DenseEquations,
        SparseEquations,
        LIS,
        PARDISO,
        PETSC
    };
};

class LinearEquationsFactory
{
public:
    static ILinearEquations* create(LinearEquationsType::type eq_type)
    {
        switch (eq_type) {
        case LinearEquationsType::DenseEquations:
            return new DenseLinearEquations();
        case LinearEquationsType::SparseEquations:
            return new SparseLinearEquation();
        case LinearEquationsType::LIS:
            return new CRSLisSolver();
        default:
            assert(false);
            return 0;
        }
    }
};
}
