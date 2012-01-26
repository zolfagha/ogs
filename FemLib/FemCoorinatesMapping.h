
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Mapping.h"
#include "FemLib/ShapeFunction.h"

namespace FemLib
{

/**
 * \brief Properties of coordinate mapping at particular locations 
 */
struct CoordMappingProperties
{
    ///// mapping locations in computed coordinates
    //double r[3];
    /// shape function N(r)
    MathLib::Matrix<double> *shape_r;
    /// gradient of shape functions, dN(r)/dx
    MathLib::Matrix<double> *dshape_dr;
    /// Jacobian matrix, J=dx/dr
    MathLib::Matrix<double> *jacobian_dxdr;
    /// determinant of the Jacobian
    double det_jacobian;
    /// inverse of the Jacobian
    MathLib::Matrix<double> *inv_jacobian_drdx;

    CoordMappingProperties()
    {
        shape_r = new MathLib::Matrix<double>(); 
        dshape_dr = new MathLib::Matrix<double>(); 
        jacobian_dxdr = new MathLib::Matrix<double>(); 
        inv_jacobian_drdx = new MathLib::Matrix<double>(); 
    }

    ~CoordMappingProperties()
    {
        delete shape_r;
        delete dshape_dr;
        delete jacobian_dxdr;
        delete inv_jacobian_drdx;
    }
};

/**
 * \brief IFemCoordinatesMapping is an interface class for geometrical mapping between actual coordinates and computational coordinates.
 */
class IFemCoordinatesMapping
{
public:
    /// initialize element
    virtual void initialize(MeshLib::IElement* ele) = 0;
    /// compute shape functions
    virtual const CoordMappingProperties* compute(const double* natural_pt) = 0;
    /// map natural coordinates to physical coordinates
    virtual double* mapToPhysicalCoordinates(const double* natural_pt) = 0;
    /// map physical coordinates to natural coordinates
    virtual double* mapFromPhysicalCoordinates(const double* physical_pt) = 0;

};

/**
 * \brief Mapping element shapes to natural coordinates
 *
 * FemNaturalCoordinates mapping converts element shapes in physical coordinates (x,y,z) to that in natural coordinates (r,s,t).
 * - Given physical coordinates should correspond to dimensions of the element, i.e (x,y) for triangles
 * - x(r,s,t) = N(r,s,t) x_i
 */
class FemNaturalCoordinates : public IFemCoordinatesMapping
{
private:
    CoordMappingProperties* _prop;
    IFemShapeFunction* _shape;
    MeshLib::IElement* _ele;

public:
    FemNaturalCoordinates() 
    {
        _prop = new CoordMappingProperties();
    }
    virtual ~FemNaturalCoordinates()
    {
        delete _prop;
    }

    /// initialize element
    virtual void initialize(MeshLib::IElement* ele)
    {
        _ele = ele;
        _shape = getShapeFunction();
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();
        _prop->shape_r->resize(1, nnodes);
        _prop->dshape_dr->resize(dim, nnodes);
        _prop->jacobian_dxdr->resize(dim, dim);
        _prop->inv_jacobian_drdx->resize(dim, dim);
    }

    /// compute mapping properties at the given location in natural coordinates
    virtual const CoordMappingProperties* compute(const double* natural_pt) 
    {
        //prepare
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();
        MathLib::Matrix<double> x(dim, nnodes);
        for (size_t i=0; i<nnodes; i++)
            for (size_t j=0; j<dim; j++)
                x(j,i) = _ele->getNodeCoordinates(i)->getData()[j]; //TODO should be via local coordinates

        //shape, dshape
        MathLib::Matrix<double> *shape  = _prop->shape_r;
        MathLib::Matrix<double> *dshape_dr = _prop->dshape_dr;
        _shape->computeShapeFunction(natural_pt, (double*)shape->getData());
        _shape->computeGradShapeFunction(natural_pt, (double*)dshape_dr->getData());

        //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
        MathLib::Matrix<double> *jac = _prop->jacobian_dxdr;
        for (size_t i_r=0; i_r<dim; i_r++) {
            for (size_t j_x=0; j_x<dim; j_x++) {
                for (size_t k=0; k<nnodes; j_x++) {
                    (*jac)(i_r,j_x) += (*dshape_dr)(i_r, k) * x(j_x,k);
                }
            }
        }

        //determinant of J
        //TODO add a method to calculate a determinant
        double det_j = 0;
        if (dim==1) {
            det_j = (*jac)(0,0);
        } else if (dim==2) {
            det_j = (*jac)(0,0)*(*jac)(1,1);
            det_j -= (*jac)(0,1)*(*jac)(1,0);
        } else if (dim==3) {
            det_j = (*jac)(0,0)*((*jac)(1,1)*(*jac)(2,2)-(*jac)(2,1)*(*jac)(1,2));
            det_j -= (*jac)(1,0)*((*jac)(0,1)*(*jac)(2,2)-(*jac)(2,1)*(*jac)(0,2));
            det_j += (*jac)(2,0)*((*jac)(0,1)*(*jac)(1,2)-(*jac)(1,1)*(*jac)(0,2));
        }
        _prop->det_jacobian = det_j;

        //inverse of J
        //TODO add a method to get inverse of the J
        MathLib::Matrix<double> *inv_jac = _prop->inv_jacobian_drdx;
        if (dim==1) {
            (*inv_jac)(0,0) = (*jac)(0,0);
        } else if (dim==2) {
            (*inv_jac)(0,0) = (*jac)(1,1);
            (*inv_jac)(0,1) = -(*jac)(0,1);
            (*inv_jac)(1,0) = -(*jac)(1,0);
            (*inv_jac)(1,1) = (*jac)(0,0);
        } else if (dim==3) {
            (*inv_jac)(0,0) = (*jac)(1,1)*(*jac)(2,2)-(*jac)(2,1)*(*jac)(1,2);
            (*inv_jac)(0,1) = (*jac)(0,2)*(*jac)(2,1)-(*jac)(0,1)*(*jac)(2,2);
            (*inv_jac)(0,2) =  (*jac)(0,1)*(*jac)(1,2)-(*jac)(0,2)*(*jac)(1,1);
            //
            (*inv_jac)(1,0) =  (*jac)(1,2)*(*jac)(2,0)-(*jac)(2,2)*(*jac)(1,0);
            (*inv_jac)(1,1) =  (*jac)(0,0)*(*jac)(2,2)-(*jac)(2,0)*(*jac)(0,2);
            (*inv_jac)(1,2) =  (*jac)(0,2)*(*jac)(1,0)-(*jac)(1,2)*(*jac)(0,0);
            //
            (*inv_jac)(2,0) =  (*jac)(1,0)*(*jac)(2,1)-(*jac)(2,0)*(*jac)(1,1);
            (*inv_jac)(2,1) =  (*jac)(0,1)*(*jac)(2,0)-(*jac)(2,1)*(*jac)(0,0);
            (*inv_jac)(2,2) =  (*jac)(0,0)*(*jac)(1,1)-(*jac)(1,0)*(*jac)(0,1);
        }
        (*inv_jac) /= det_j;

        return _prop;
    };

    double* mapToPhysicalCoordinates(const double* natural_pt);
    double* mapFromPhysicalCoordinates(const double* physical_pt);

private:
    virtual IFemShapeFunction* getShapeFunction() = 0;
};

/**
 * Axisymmetric class
 * - coordinates (r, theta, z)
 * - dxdydz = 2pi*r*drdz
 */
class FemAxisymmetric
{
public:
    FemAxisymmetric() : _r(.0) {};
    virtual ~FemAxisymmetric() {};

    virtual void computeMappingFunctions(double* natural_pt, FemMapComputation::type compType) 
    {
        //FemIsoparametricMapping::computeMappingFunctions(natural_pt, compType);
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
        //MeshLib::IElement *e = getElement();
        //double *shapefct = getShapeFunction();
        //for (size_t i=0; i<e->getNumberOfNodes(); i++) {
        //    const GeoLib::Point *pt = e->getMappedGeometry()->getNodePoint(i);
        //    r += shapefct[i] * (*pt)[0];
        //}
        return r;
    }
};

class FemLowerDimension
{

};

}
