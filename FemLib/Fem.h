
#pragma once

#include "MeshLib/Core/IElement.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "Integration.h"

namespace FemLib
{

struct FiniteElementType
{
    enum type {
        Lagrange = 1,
        Hermite = 2,
        INVALID = -1
    };
};


class IFiniteElement
{
public:
    //virtual void configure(MeshLib::IElement* ele) = 0;

    //virtual void assembleNN(MathLib::Matrix _m, double func) = 0;
    //virtual void assembledNdN(MathLib::Matrix _m, double func) = 0;
};


class FemLagrangeElement : public IFiniteElement
{
public:
    FemLagrangeElement() {};

    virtual void configure(MeshLib::IElement* ele, int degree) {
        // set mapping
        // calculate shape functions
        // 
    }

    IFemIntegration* getIntegration() {return &_integration;};
    void getShapeFunctions(int igp, double* shape, double* dshape, double *jacobian);
    double* getShapeFunction(int igp);
    MathLib::Matrix<double>* getGradShapeFunction(int igp);
    double* getTestFunction(int igp);
    MathLib::Matrix<double>* getGradTestFunction(int igp);
    double getDetJacobian(int igp);

    double* computeShapeFunction( double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }
    double* computeTestFunction( double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    MathLib::Matrix<double>* computeGradShapeFunction( double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }
    MathLib::Matrix<double>* computeGradTestFunction( double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    void integrate( void (*calcLap)(double*, MathLib::Matrix<double>&), MathLib::Matrix<double> M ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }




private:
    IFemIntegration _integration;

    int _ele_geo;
    int _shape_func;
};

class FiniteElementFactory
{
public:
    IFiniteElement* createFiniteElement(FiniteElementType::type ele_type) {
        switch (ele_type) {
        case FiniteElementType::Lagrange:
            return new FemLagrangeElement();
            break;
        case FiniteElementType::Hermite:
            break;
        default:
            break;
        }

        return 0;
    }    
};

}
