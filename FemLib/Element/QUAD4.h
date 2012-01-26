#pragma once

#include "IFemElement.h"

namespace FemLib
{

class QUAD4 : public FeBaseNaturalCoordinates<FiniteElementType::QUAD4, 4>
{
private:
    virtual IFemCoordinatesMapping* getMappingMethod() const 
    {
        return new FemNaturalCoordinates();
    }
    virtual IFemShapeFunction* getBaseFunction() const 
    {
        return new FemShapeQuad4();
    }

    virtual IFemIntegration* getIntegrationMethod() const
    {
        return new FemIntegrationGauss();
    }


    ///// initialize object for given mesh elements
    //virtual void configure( const MeshLib::IElement * e )
    //{

    //}
    ///// make interpolation from nodal values
    //virtual double interpolate(double *pt, double *nodal_values) 
    //{
    //    return .0;
    //}
    ///// extrapolate nodal values from integration point values
    //virtual void extrapolate(double *integral_values, double *nodal_values) 
    //{

    //}
    ///// compute an matrix M = Int{W^T F N} dV
    //virtual void computeIntTestShape( void (*func)(double*), MathLib::Matrix<double> *) 
    //{

    //}
    ///// compute an matrix M = Int{W^T F dN} dV
    //virtual void computeIntTestDShape( void (*func)(double*), MathLib::Matrix<double> *)
    //{

    //}
    ///// compute an matrix M = Int{dW^T F dN} dV
    //virtual void computeIntDTestDShape( void (*func)(double*), MathLib::Matrix<double> *)
    //{

    //}
    ///// compute an vector V = Int{W^T*N}dV*U_i
    //virtual void computeIntTestShapeNodalVal(double *nod_val, double *result)
    //{

    //}

    ///// return finite element type
    //virtual const FiniteElementType::type& getFeType() const 
    //{
    //    return FiniteElementType::QUAD4;
    //}

    ///// return integration method
    //virtual IFemIntegration* getIntegrationMethod() const 
    //{
    //    return 0;
    //}

    ///// return mapping method
    //virtual IFemMapping* getMapping() const
    //{
    //    return 0;
    //}

    ///// return the number of DoFs
    //virtual size_t getNumberOfDOFs() const
    //{
    //    return 4;
    //};
};

}
