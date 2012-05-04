
#pragma once

#include <vector>
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "FemDirichletBC.h"

namespace FemLib
{

class IDirichletBCMethod
{
public:
//    virtual void apply(int linearEqs, DirichletBC &bc) = 0;
};

class DiagonalizeMethod : public IDirichletBCMethod
{
public:
    void apply( MathLib::ILinearEquations& eqs, FemDirichletBC &bc)
    {
        eqs.setKnownX(bc.getListOfBCNodes(), bc.getListOfBCValues());
    }
};

class ResizeMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class PenaltyMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class LagrangeMultiplier : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};


} //end
