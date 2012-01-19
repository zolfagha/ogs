
#pragma once

#include "MeshLib/Core/IElement.h"
#include "MeshLib/Core/Element.h"
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
    virtual void configure( MeshLib::IElement * e ) = 0;
    virtual void integrate( void (*func)(double*, MathLib::Matrix<double>*), MathLib::Matrix<double>* M ) = 0;
    virtual void integrateShapeShape( void (*func)(double*), MathLib::Matrix<double> *) = 0;
    virtual void integrateShapeDShape( void (*func)(double*), MathLib::Matrix<double> *) = 0;
    virtual void integrateDShapeDShape( void (*func)(double*), MathLib::Matrix<double> *) = 0;
    virtual const MathLib::Matrix<double> * computeGradShapeFunction( const double * pt) = 0; 

    IFemIntegration * getIntegrationMethod() 
    {
        return _integration_method;
    }


    //virtual void assembleNN(MathLib::Matrix _m, double func) = 0;
    //virtual void assembledNdN(MathLib::Matrix _m, double func) = 0;
private:
    IFemIntegration* _integration_method;
};


/**
 * \ingroup FemLib
 *
 * \brief Lagrange type Finite Element
 */
class FemLagrangeElement : public IFiniteElement
{
public:
    FemLagrangeElement() {};

    virtual void configure(MeshLib::IElement* ele) {
        // set mapping
        // calculate shape functions
        // 
    }

    IFemIntegration* getIntegration() {return _integration;};
    void getShapeFunctions(int igp, double* shape, double* dshape, double *jacobian);
    double* getShapeFunction(int igp);
    MathLib::Matrix<double>* getGradShapeFunction(int igp);
    double* getTestFunction(int igp);
    MathLib::Matrix<double>* getGradTestFunction(int igp);
    double getDetJacobian(int igp);

    double* computeShapeFunction( const double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }
    double* computeTestFunction( const double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    const MathLib::Matrix<double> * computeGradShapeFunction( const double * pt)
    {
        throw std::exception("The method or operation is not implemented.");
    }

    MathLib::Matrix<double>* computeGradTestFunction( const double * pt ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    virtual void integrate( void (*func)(double*, MathLib::Matrix<double>*), MathLib::Matrix<double>* M ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    virtual void integrateShapeShape( void (*func)(double*), MathLib::Matrix<double> *) {

    }
    virtual void integrateShapeDShape( void (*func)(double*), MathLib::Matrix<double> *) {

    }
    virtual void integrateDShapeDShape( void (*func)(double*), MathLib::Matrix<double> *) {

    }



private:
    IFemIntegration *_integration;

    int _ele_geo;
    int _shape_func;
};


class FEConstantLinearTriangle
{
public:
    void assembleNN(MathLib::Matrix<double> &mat) {
        MeshLib::Triangle *ele;
        const size_t n_ele_nodes = 3;
        double a[n_ele_nodes], b[n_ele_nodes], c[n_ele_nodes];
        double nodes_x[n_ele_nodes], nodes_y[n_ele_nodes], nodes_z[n_ele_nodes];
        // xyz
        //for (size_t i=0; i<n_ele_nodes; i++) {
        //    msh->getNodeCoordinates(ele->getNodeID(i), pt);
        //    nodes_x[i] = pt[0];
        //    nodes_y[i] = pt[1];
        //    nodes_z[i] = pt[2];
        //}
        // area
        const double A = 0.5*(nodes_x[0]*(nodes_y[1]-nodes_y[2])+nodes_x[1]*(nodes_y[2]-nodes_y[0])+nodes_x[2]*(nodes_y[0]-nodes_y[1]));
        // set a,b,c
        a[0] = 0.5/A*(nodes_x[1]*nodes_y[2]-nodes_x[2]*nodes_y[1]);
        b[0] = 0.5/A*(nodes_y[1]-nodes_y[2]);
        c[0] = 0.5/A*(nodes_x[2]-nodes_x[1]);
        a[1] = 0.5/A*(nodes_x[2]*nodes_y[0]-nodes_x[0]*nodes_y[2]);
        b[1] = 0.5/A*(nodes_y[2]-nodes_y[0]);
        c[1] = 0.5/A*(nodes_x[0]-nodes_x[2]);
        a[2] = 0.5/A*(nodes_x[0]*nodes_y[1]-nodes_x[1]*nodes_y[0]);
        b[2] = 0.5/A*(nodes_y[0]-nodes_y[1]);
        c[2] = 0.5/A*(nodes_x[1]-nodes_x[0]);

        // assemble local EQS
        // Int{w S ph/pt + div(w) K div(p)}dA = Int{w K div(p)}dL
        mat(0,0) = b[0]*b[0] + c[0]*c[0];
        mat(0,1) = b[0]*b[1] + c[0]*c[1];
        mat(0,2) = b[0]*b[2] + c[0]*c[2];
        mat(1,1) = b[1]*b[1] + c[1]*c[1];
        mat(1,2) = b[1]*b[2] + c[1]*c[2];
        mat(2,2) = b[2]*b[2] + c[2]*c[2];
        //local_K *= A;
        // symmetric
        for (size_t i=0; i<n_ele_nodes; i++)
            for (size_t j=0; j<i; j++)
                mat(i,j) = mat(j,i);

    };

    void assembledNdN(MathLib::Matrix<double> &local_K) {
    };
};


//class FiniteElementFactory
//{
//public:
//    IFiniteElement* createFiniteElement(FiniteElementType::type ele_type) {
//        switch (ele_type) {
//        case FiniteElementType::Lagrange:
//            return new FemLagrangeElement();
//            break;
//        case FiniteElementType::Hermite:
//            break;
//        default:
//            break;
//        }
//
//        return 0;
//    }    
//};


/* \ingroup FemLib
    *
    * \brief Hermite type Finite Element. Hermite type uses additional degree 
    * of freedoms to consider continuity of derivatives.
    */
class FemHermiteElement : public IFiniteElement
{
};

}
