
#pragma once

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/Element.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "Integration.h"
#include "FemCoorinatesMapping.h"
#include "ShapeFunction.h"

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

typedef double (*Fscalar)(const double *x);
typedef double* (*Fvector)(const double *x);

/**
 * \brief IFiniteElement class is an interface to all kinds of finite element classes. 
 *
 * Roles of the class is
 * - assembly of typical local matrices, e.g. N^T N, dN^T dN
 * - interpolation, extrapolation
 * - integration over domain, boundary
 */
class IFiniteElement
{
public:
    /// initialize object for given mesh elements
    virtual void configure( MeshLib::IMesh * msh, MeshLib::IElement * e ) = 0;
    /// return finite element type
    virtual const FiniteElementType::type getFeType() const = 0;
    /// return mesh element
    virtual MeshLib::IElement* getElement() const = 0;
    /// return the number of variables
    virtual size_t getNumberOfVariables() const = 0;

    /// compute basis functions at the given point
    virtual void computeBasisFunctions(const double *x) = 0;
    /// get evaluated basis functions. 
    virtual MathLib::Matrix<double>* getBasisFunction() = 0;
    /// get evaluated gradient of basis functions. 
    virtual MathLib::Matrix<double>* getGradBasisFunction() = 0;

    /// make interpolation from nodal values
    virtual double interpolate(double *pt, double *nodal_values) = 0;
    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN( Fscalar, MathLib::Matrix<double> &) = 0;
    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN( Fvector, MathLib::Matrix<double> &) = 0;
    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN( Fscalar, MathLib::Matrix<double> &) = 0;

    /// get the integration method
    virtual IFemNumericalIntegration* getIntegrationMethod() const = 0;
};

/**
 * \brief Base class for implementation of finite element classes
 */
template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES>
class TemplateFeBase : public IFiniteElement
{
public:
    TemplateFeBase() : _ele(0) {};
    virtual ~TemplateFeBase() {};

    /// return mesh element
    MeshLib::IElement* getElement() const {return _ele;};
    /// return the number of variables
    size_t getNumberOfVariables() const { return N_VARIABLES; };
    /// return finite element type
    const FiniteElementType::type getFeType() const { return T_FETYPE; };
protected:
    void setElement(MeshLib::IElement* e) {_ele = e;};
private:
    MeshLib::IElement* _ele;
};


/**
 * \brief Base for any isoparametric FE classes
 */
template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES, class T_SHAPE, class T_INTEGRAL>
class FeBaseIsoparametric : public TemplateFeBase<T_FETYPE, N_VARIABLES>
{
public:
    FeBaseIsoparametric()
    {
        _mapping = new FemNaturalCoordinates(new T_SHAPE());
        _integration = new T_INTEGRAL();
    };

    virtual ~FeBaseIsoparametric() 
    {
        delete _mapping;
        delete _integration;
    };

    virtual IFemNumericalIntegration* getIntegrationMethod() const {return _integration;};

    /// initialize object for given mesh elements
    virtual void configure( MeshLib::IMesh * msh, MeshLib::IElement * e )
    {
        setElement(e);
        _mapping->initialize(e);
    }

    virtual void computeBasisFunctions(const double *x)
    {
        _mapping->compute(x);
    }

    virtual MathLib::Matrix<double>* getBasisFunction()
    {
        return _mapping->getProperties()->shape_r;
    }

    virtual MathLib::Matrix<double>* getGradBasisFunction()
    {
        return _mapping->getProperties()->dshape_dx;
    }


    /// make interpolation from nodal values
    virtual double interpolate(double *natural_pt, double *nodal_values)
    {
        const CoordMappingProperties *prop = _mapping->compute(natural_pt);
        double *N = (double*)prop->shape_r->getData();
        double v = .0;
        for (size_t i=0; i<getNumberOfVariables(); i++)
            v+=N[i]*nodal_values[i];
        return v;
    }

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(Fscalar f, MathLib::Matrix<double> &mat)
    {
        const size_t n_gp = _integration->getNumberOfSamplingPoints();
        double x[3];
        for (size_t i=0; i<n_gp; i++) {
            _integration->getSamplingPoint(i, x);
            const CoordMappingProperties *coord_prop = _mapping->compute(x);
            MathLib::Matrix<double> *basis = coord_prop->shape_r;
            MathLib::Matrix<double> *test = coord_prop->shape_r;
            double fac = coord_prop->det_jacobian * _integration->getWeight(i);
            double v = f(x);
            fac *= v;
            test->transposeAndMultiply(*basis, mat, fac);
        }
    }

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(Fvector f, MathLib::Matrix<double> &mat)
    {
        const size_t n_gp = _integration->getNumberOfSamplingPoints();
        double x[3];
        for (size_t i=0; i<n_gp; i++) {
            _integration->getSamplingPoint(i, x);
            const CoordMappingProperties *coord_prop = _mapping->compute(x);
            MathLib::Matrix<double> *dbasis = coord_prop->dshape_dx;
            MathLib::Matrix<double> *test = coord_prop->dshape_dx;
            double fac = coord_prop->det_jacobian * _integration->getWeight(i);
            double *v = f(x);
            test->transposeAndMultiply(*dbasis, v, mat, fac);
        }
    }

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(Fscalar f, MathLib::Matrix<double> &mat)
    {
        const size_t n_gp = _integration->getNumberOfSamplingPoints();
        double x[3];
        for (size_t i=0; i<n_gp; i++) {
            _integration->getSamplingPoint(i, x);
            const CoordMappingProperties *coord_prop = _mapping->compute(x);
            MathLib::Matrix<double> *dbasis = coord_prop->dshape_dx;
            MathLib::Matrix<double> *dtest = coord_prop->dshape_dx;
            double fac = coord_prop->det_jacobian * _integration->getWeight(i);
            double v = f(x);
            fac *= v;
            dtest->transposeAndMultiply(*dbasis, mat, fac);
        }
    }

private:
    FemNaturalCoordinates *_mapping;
    IFemNumericalIntegration* _integration;

//protected:
//    virtual IFemShapeFunction* createShapeFunction() const = 0;
//    virtual IFemIntegration* createIntegrationMethod() const = 0;
};


}

