/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LinearEquationsFactory.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

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
