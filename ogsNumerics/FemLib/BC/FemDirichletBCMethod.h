
#pragma once

#include <vector>
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

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
    void apply(const std::vector<size_t> &bc_nodes, const std::vector<double> &bc_values, MathLib::ILinearEquations& eqs)
    {
        eqs.setKnownX(bc_nodes, bc_values);
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
