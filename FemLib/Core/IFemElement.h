
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/Element.h"

#include "Integration.h"
#include "FemCoorinatesMapping.h"
#include "ShapeFunction.h"
#include "FemLib/Core/Element/FemElementList.h"

namespace FemLib
{

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
    /// setup object for given mesh elements
    virtual void configure(MeshLib::IElement * e ) = 0;
    /// return finite element type
    virtual const FiniteElementType::type getFeType() const = 0;
    /// return this mesh element
    virtual MeshLib::IElement* getElement() const = 0;
    /// return the number of variables
    virtual size_t getNumberOfVariables() const = 0;

    /// compute basis functions \f$ \mathbf{N}_e \f$, \f$ {\nabla}_x \mathbf{N}_e \f$ at the given point
    virtual void computeBasisFunctions(const double *x) = 0;
    /// get evaluated basis functions \f$ \mathbf{N}_e \f$. 
    virtual MathLib::Matrix<double>* getBasisFunction() = 0;
    /// get evaluated gradient of basis functions \f$ {\nabla}_x \mathbf{N}_e \f$. 
    virtual MathLib::Matrix<double>* getGradBasisFunction() = 0;

    /// make interpolation from nodal values \f$ u^h(\mathbf{x}) = \mathbf{N}_e (\mathbf x) \mathbf{u}_e \f$
    virtual double interpolate(double *pt, double *nodal_values) = 0;
    /// compute an matrix \f$ \mathbf{M}_e = \int_{\Omega_e} {\mathbf{N}_e^*}^T f(\mathbf x) \mathbf{N}_e d\Omega \f$
    virtual void integrateWxN( MathLib::IFunction<double, double*>*, MathLib::Matrix<double> &) = 0;
    /// compute an matrix \f$ \mathbf{M}_e = \int_{\Omega_e} {\mathbf{N}_e^*}^T \mathbf{f}(\mathbf x) \nabla \mathbf{N}_e d\Omega \f$
    virtual void integrateWxDN( MathLib::IFunction<double*, double*>*, MathLib::Matrix<double> &) = 0;
    /// compute an matrix \f$ \mathbf{M}_e = \int_{\Omega_e} {\nabla \mathbf{N}_e^*}^T f(\mathbf x) \nabla \mathbf{N}_e d\Omega \f$
    virtual void integrateDWxDN( MathLib::IFunction<double, double*> *f, MathLib::Matrix<double> &) = 0;

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
    TemplateFeBase(MeshLib::IMesh *msh) : _msh(msh), _ele(0) {};
    virtual ~TemplateFeBase() {};

    void setMesh(MeshLib::IMesh * msh) {_msh = msh;};
    const MeshLib::IMesh* getMesh() const {return _msh;};
    /// return mesh element
    MeshLib::IElement* getElement() const {return _ele;};
    /// return the number of variables
    size_t getNumberOfVariables() const { return N_VARIABLES; };
    /// return finite element type
    const FiniteElementType::type getFeType() const { return T_FETYPE; };
protected:
    void setElement(MeshLib::IElement* e) {_ele = e;};
private:
    MeshLib::IMesh* _msh;
    MeshLib::IElement* _ele;
};


}

