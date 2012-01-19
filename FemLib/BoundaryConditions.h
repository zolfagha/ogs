
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
class DirichletBC //per variable?
{
public:
    template<typename Tval, typename Tpos>
    void set(FEMNodalFunction<Tval,Tpos> *fem, GeoLib::GeoObject *geo, MathLib::IDistribution *func)
    {
        throw std::exception("The method or operation is not implemented.");
    }

    template<typename Tmat, typename Tvec>
    void apply( Tmat* globalA, Tvec* globalRHS ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }


private:
    // node id, var id, value
    //IDirichletBCMethod *_method;
};

class NeumannBC
{
public:
    template<typename Tval, typename Tpos>
    void set(FEMNodalFunction<Tval,Tpos> *fem, GeoLib::GeoObject *geo, MathLib::IDistribution *func) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    size_t getNumberOfConditions() const
    {
        throw std::exception("The method or operation is not implemented.");
    }

    size_t getConditionDoF( size_t i ) const
    {
        throw std::exception("The method or operation is not implemented.");
    }

    double getConditionValue( size_t i ) const
    {
        throw std::exception("The method or operation is not implemented.");
    }

    template<typename T>
    void apply( T* globalRHS ) 
    {
        for (size_t i=0; i<this->getNumberOfConditions(); i++)
            (*globalRHS)[this->getConditionDoF(i)] += this->getConditionValue(i);
    }



private:
    // node id, var id, value
};


class CauchyBC
{
public:
    void set(FEMNodalFunction<double,double> *fem, int geo, int func);
};

class OpenBC
{
public:
    void set(FEMNodalFunction<double,double> *fem, int geo, int func);
};


//----------------------------------------------------------

class IDirichletBCMethod
{
public:
//    virtual void apply(int linearEqs, DirichletBC &bc) = 0;
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
