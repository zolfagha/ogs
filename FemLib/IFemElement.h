
#pragma once

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/Element.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "Integration.h"
#include "FemCoorinatesMapping.h"
#include "ShapeFunction.h"
#include "TestFunction.h"

namespace FemLib
{

/// Finite element type
struct FiniteElementType
{
    enum type {
        LINE2,
        LINE3,
        TRI3,
        TRI3CONST,
        TRI6,
        QUAD4,
        QUAD8,
        QUAD9,
        INVALID
    };
};

typedef double (*Fscalar)(double *x);
typedef double* (*Fvector)(double *x);

/**
 * IFiniteElement class is an interface to all kinds of finite element classes. 
 * Roles of the class is
 * - assembly of typical local matrices, e.g. N^T N, dN^T dN
 * - interpolation, extrapolation
 * - integration over domain, boundary
 */
class IFiniteElement
{
public:
    IFiniteElement() : _ele(0) {};
    virtual ~IFiniteElement() {};

    /// initialize object for given mesh elements
    virtual void configure( MeshLib::IMesh * msh, MeshLib::IElement * e ) = 0;
    /// return finite element type
    virtual const FiniteElementType::type getFeType() const = 0;
    /// return mesh element
    MeshLib::IElement* getElement() const {return _ele;};
    /// return the number of variables
    virtual size_t getNumberOfVariables() const = 0;

    /// make interpolation from nodal values
    virtual double interpolate(double *pt, double *nodal_values) = 0;
    /// compute an matrix M = Int{W^T F N} dV
    virtual void computeIntTestShape( Fscalar, MathLib::Matrix<double> &) = 0;
    /// compute an matrix M = Int{W^T F dN} dV
    virtual void computeIntTestDShape( Fvector, MathLib::Matrix<double> &) = 0;
    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void computeIntDTestDShape( Fscalar, MathLib::Matrix<double> &) = 0;

private:
    MeshLib::IElement* _ele;
};

template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES>
class TemplateFeBase : public IFiniteElement
{
public:
    /// return the number of variables
    size_t getNumberOfVariables() const { return N_VARIABLES; };
    /// return finite element type
    const FiniteElementType::type getFeType() const { return T_FETYPE; };
};

// use of natural coordinates

template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES>
class FeBaseNaturalCoordinates : public TemplateFeBase<T_FETYPE, N_VARIABLES>
{
public:
    FeBaseNaturalCoordinates() : _ele(0) 
    {
        _base = getBaseFunction();
        _test = getTestFunction();
        _integral = getIntegrationMethod();
        _mapping = new FemNaturalCoordinates(_base);
    };

    virtual ~FeBaseNaturalCoordinates() 
    {
        delete _mapping;
        delete _base;
        delete _test;
        delete _integral;
    };

    /// initialize object for given mesh elements
    virtual void configure( const MeshLib::IElement * e )
    {

    }


    /// make interpolation from nodal values
    virtual double interpolate(double *pt, double *nodal_values)
    {
        double v = 0;
        MathLib::Matrix<double> *matN = _base->computeShapeFunction(pt);
        double *N = (double*)matN->getData();
        for (size_t i=0; i<getNumberOfVariables(); i++)
            v+=N[i]*nodal_values[i];
        return v;
    }

    /// extrapolate nodal values from integration point values
    virtual void extrapolate(double *integral_values, double *nodal_values)
    {
        throw std::exception("not implemented yet.");
    }

    /// compute an matrix M = Int{W^T F N} dV
    virtual void computeIntTestShape( void (*func)(double*), MathLib::Matrix<double> *mat)
    {
    }

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void computeIntTestDShape( void (*func)(double*), MathLib::Matrix<double> *)
    {

    }

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void computeIntDTestDShape( void (*func)(double*), MathLib::Matrix<double> *)
    {

    }

    /// compute an vector V = Int{W^T*N}dV*U_i
    virtual void computeIntTestShapeNodalVal(double *nod_val, double *result) 
    {
        MathLib::Matrix<double> M(this->getNumberOfVariables(), this->getNumberOfVariables());
        //for each sampling point
       size_t n_gp = _integral->getNumberOfSamplingPoints();
       for (size_t i=0; i<n_gp; i++) {
            const double *x = _integral->getSamplingPoint(i);
            //TODO Mapping should provide N, dN, j because dr/dx is involved
            MathLib::Matrix<double> *N = _base->computeShapeFunction(x);
            MathLib::Matrix<double> *T = _test->computeTestFunction(*N);
            _mapping->compute(x, FemMapComputation::ALL);
            double fac = _mapping->getDetJacobian() * _integral->getWeight(i);
            //M += W^T N * j *w
            T->transposeAndMultiply(*N, M, fac);
       }
       //r = M * u
       M.axpy(1.0, nod_val, 0.0, result);
    }

private:
    MeshLib::IElement* _ele;
    IFemCoordinatesMapping *_mapping;
    IFemShapeFunction *_base;
    IFemTestFunction *_test;
    IFemIntegration *_integral;

    virtual IFemShapeFunction* getBaseFunction() const = 0;
    virtual IFemIntegration* getIntegrationMethod() const = 0;

    virtual IFemTestFunction* getTestFunction() const {
        return new FemTestFunctionBubnov();
    }
};

}

