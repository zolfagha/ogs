
#pragma once

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/Element.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "Integration.h"
#include "Mapping.h"

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
    virtual void configure( const MeshLib::IElement * e ) = 0;
    /// make interpolation from nodal values
    virtual double interpolate(double *pt, double *nodal_values) = 0;
    /// extrapolate nodal values from integration point values
    virtual void extrapolate(double *integral_values, double *nodal_values) = 0;
    /// compute an matrix M = Int{W^T F N} dV
    virtual void computeIntTestShape( void (*func)(double*), MathLib::Matrix<double> *) = 0;
    /// compute an matrix M = Int{W^T F dN} dV
    virtual void computeIntTestDShape( void (*func)(double*), MathLib::Matrix<double> *) = 0;
    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void computeIntDTestDShape( void (*func)(double*), MathLib::Matrix<double> *) = 0;
    /// compute an vector V = Int{W^T*N}dV*U_i
    virtual void computeIntTestShapeNodalVal(double *nod_val, double *result) = 0;

    /// return finite element type
    virtual const FiniteElementType::type& getFeType() const = 0;
    /// return mesh element
    MeshLib::IElement* getElement() const {return _ele;};
    /// return integration method
    virtual IFemIntegration* getIntegrationMethod() const = 0;
    /// return mapping method
    virtual IFemMapping* getMapping();
    /// return the number of DoFs
    virtual size_t getNumberOfDOFs() const = 0;

private:
    MeshLib::IElement* _ele;
};

}
