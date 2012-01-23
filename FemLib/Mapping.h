
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Mapping.h"

namespace FemLib
{

/// 
struct FemMapComputation
{
    enum type {
        SHAPE,
        DSHAPE,
        ALL
    };
};

/**
 * IFemMapping is an interface class for geometrical mapping between actual coordinates and computational coordinates.
 */
class IFemMapping
{
public:
    /// initialize element
    virtual void configure(const MeshLib::IElement* ele) = 0;
    /// compute shape functions
    virtual void computeMappingFunctions(const double* natural_pt, FemMapComputation::type compType) = 0;
    /// return shape functions
    virtual double* getShapeFunction() = 0;
    /// return derivatives of shape functions
    virtual const MathLib::Matrix<double>* getGradShapeFunction() = 0;
    /// map natural coordinates to physical coordinates
    virtual double* mapNatural2Physical(const double* natural_pt) = 0;
    /// map physical coordinates to natural coordinates
    virtual double* mapPhysical2Natural(const double* physical_pt) = 0;
};

/**
 * 
 */
class FemInvariantMapping : public IFemMapping
{
public:
};

/**
 * Isoparametric mapping converts element shapes in physical coordinates (x,y,z) to that in natural coordinates (r,s,t).
 * - Given physical coordinates should correspond to dimensions of the element, i.e (x,y) for triangles
 * - x(r,s,t) = N(r,s,t) x_i
 */
class FemIsoparametricMapping : public IFemMapping
{
public:
    ///
    virtual void configure(MeshLib::IElement* ele) 
    {
        _ele = ele;
    };

    /// 
    virtual void computeMappingFunctions(double* natural_pt, FemMapComputation::type compType) 
    {
        if (compType==FemMapComputation::SHAPE) {
            computeShapeFunction(natural_pt);
        } else if (compType==FemMapComputation::DSHAPE) {
            computeGradShapeFunction(natural_pt);
        } else if (compType==FemMapComputation::ALL) {
            computeShapeFunction(natural_pt);
            computeGradShapeFunction(natural_pt);
        }
    };

    ///
    MeshLib::IElement* getElement() const {return _ele;};
    double* mapNatural2Physical(double* natural_pt);
    double* mapPhysical2Natural(double* physical_pt);

    double* getShapeFunction();
    void getJacobian();
    void getInvJacobian();
    double getDetJacobian();
private:
    void computeShapeFunction(double*);
    void computeGradShapeFunction(double*);
    void computeJacobian(double*);

    MeshLib::IElement* _ele;
    double _det_jac;
    double *_jac;
    double *_inv_jac;
    double *_localX;
    double *_localY;
    double *_localZ;
};


/**
 * Axisymmetric class
 * - coordinates (r, theta, z)
 * - dxdydz = 2pi*r*drdz
 */
class FemAxisymmetric : public FemIsoparametricMapping
{
public:
    FemAxisymmetric() : _r(.0) {};
    virtual ~FemAxisymmetric() {};

    virtual void computeMappingFunctions(double* natural_pt, FemMapComputation::type compType) 
    {
        FemIsoparametricMapping::computeMappingFunctions(natural_pt, compType);
        _r = this->computeRadius();
    };


    /// return radius
    double getRadius() const {
        return _r;
    }
    
private:
    double _r;

    double computeRadius() {
        double r = .0;
        MeshLib::IElement *e = getElement();
        double *shapefct = getShapeFunction();
        for (size_t i=0; i<e->getNumberOfNodes(); i++) {
            const GeoLib::Point *pt = e->getMappedGeometry()->getNodePoint(i);
            r += shapefct[i] * (*pt)[0];
        }
        return r;
    }
};

class FemLowerDimension
{

};

}
