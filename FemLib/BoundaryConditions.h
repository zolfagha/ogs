
#pragma once

#include "FemFunction.h"

#include "MathLib/LinearInterpolation.h"

namespace MathLib
{

class IDistribution
{
public:
    virtual double getValue(double* x) = 0;
};

class DistributionConstant : public IDistribution
{
public:
    DistributionConstant(double v) {_v = v;};
    virtual double getValue(double* x) {return _v;}
private:
    double _v;
};

class DistributionLineaer : public IDistribution
{
public:
    DistributionLineaer(double v) {};
    virtual double getValue(double* x) {return 0;}
private:
    LinearInterpolation *_linear;
};

}


namespace FemLib
{

//----------------------------------------------------------
class CauchyBC
{
public:
    void set(FEMNodalFunction<double,double> *fem, int geo, int func);
};


//----------------------------------------------------------
class NeumannBC
{
public:
    void set(FEMNodalFunction<double,double> *fem, int geo, int func);
private:
    // node id, var id, value
};


//----------------------------------------------------------

class IDirichletBCMethod
{
public:
//    virtual void apply(int linearEqs, DirichletBC &bc) = 0;
};

class DirichletBC //per variable?
{
public:
    template<typename Tval, typename Tpos>
    void set(FEMNodalFunction<Tval,Tpos> *fem, GeoLib::GeoObject *geo, MathLib::IDistribution *func);
private:
    // node id, var id, value
    IDirichletBCMethod *_method;
};


class DiagonizeMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, DirichletBC &bc);
};

class ResizeMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, DirichletBC &bc);
};

class PenaltyMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, DirichletBC &bc);
};

class LagrangeMultiplier : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, DirichletBC &bc);
};

}
