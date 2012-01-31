
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "FemLib/Core/IFemElement.h"

namespace FemLib
{

/**
 * \brief TRI3CONST finite element class
 *
 * TRI3CONST is a finite element class for specific case of linear triangle elements with
 * element constant parameters. Efficient integration can be expected due to use of analytical
 * integration methods.
 * 
 * Remark:
 * - Given coefficients are assumed to be constant over a element.
 * - Axisymmetric model is not supported!
 */    
class TRI3CONST : public TemplateFeBase<FiniteElementType::TRI3CONST, 3>
{
private:
    MeshLib::Triangle *_ele;
    MeshLib::IMesh *_msh;
    double a[3], b[3], c[3];
    double A;
    FemIntegrationAnalytical _integration;
    MathLib::Matrix<double> _shape, _dshape;

    void computeBasisFunction(const double *x,  double *shape);
    void computeGradBasisFunction(const double *x,  MathLib::Matrix<double> &mat);
public:
    TRI3CONST() : _shape(1,3), _dshape(2,3) {};

    /// initialize object for given mesh elements
    void configure( MeshLib::IMesh * msh, MeshLib::IElement * e );

    /// 
    void computeBasisFunctions(const double *x);
    MathLib::Matrix<double>* getBasisFunction();
    MathLib::Matrix<double>* getGradBasisFunction();


    /// make interpolation from nodal values
    double interpolate(double *x, double *nodal_values);

    /// compute an matrix M = Int{W^T F N} dV
    void integrateWxN( MathLib::IFunction<double, double*>* f, MathLib::Matrix<double> &mat);

    /// compute an matrix M = Int{W^T F dN} dV
    void integrateWxDN( MathLib::IFunction<double*, double*>* f, MathLib::Matrix<double> &mat);

    /// compute an matrix M = Int{dW^T F dN} dV
    void integrateDWxDN( MathLib::IFunction<double, double*> *f, MathLib::Matrix<double> &mat);

    /// get the integration method
    IFemNumericalIntegration* getIntegrationMethod() const {return (IFemNumericalIntegration*)&_integration;};

};

}
