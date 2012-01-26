
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "IFemElement.h"

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

public:
    /// initialize object for given mesh elements
    void configure( MeshLib::IMesh * msh, MeshLib::Triangle * e );

    /// make interpolation from nodal values
    double interpolate(double *x, double *nodal_values);

    /// compute an matrix M = Int{W^T F N} dV
    void computeIntTestShape( Fscalar f, MathLib::Matrix<double> &mat);

    /// compute an matrix M = Int{W^T F dN} dV
    void computeIntTestDShape( Fvector f, MathLib::Matrix<double> &mat);

    /// compute an matrix M = Int{dW^T F dN} dV
    void computeIntDTestDShape( Fscalar f, MathLib::Matrix<double> &mat);
};

}
