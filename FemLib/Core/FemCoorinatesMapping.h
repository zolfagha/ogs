
#pragma once

#include <iostream>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Mapping.h"
#include "FemLib/Core/ShapeFunction.h"

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
    /// gradient of shape functions, dN(r)/dr
    MathLib::Matrix<double> *dshape_dr;
    /// gradient of shape functions, dN(r)/dx
    MathLib::Matrix<double> *dshape_dx;
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
        dshape_dx = new MathLib::Matrix<double>(); 
        jacobian_dxdr = new MathLib::Matrix<double>(); 
        inv_jacobian_drdx = new MathLib::Matrix<double>(); 
    }

    ~CoordMappingProperties()
    {
        delete shape_r;
        delete dshape_dr;
        delete dshape_dx;
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
    virtual void mapToPhysicalCoordinates(const double* natural_pt, double*) = 0;
    /// map physical coordinates to natural coordinates
    virtual void mapFromPhysicalCoordinates(const double* physical_pt, double*) = 0;

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
    MathLib::Matrix<double> x;

public:
    FemNaturalCoordinates(IFemShapeFunction *shape) 
    {
        _shape = shape;
        _prop = new CoordMappingProperties();
    }
    virtual ~FemNaturalCoordinates()
    {
        delete _prop;
    }

    /// initialize element
    virtual void initialize(MeshLib::IElement* ele)
    {
        assert(ele->getMappedCoordinates()!=0);

        _ele = ele;
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();
        _prop->shape_r->resize(1, nnodes);
        _prop->dshape_dr->resize(dim, nnodes);
        _prop->dshape_dx->resize(dim, nnodes);
        _prop->jacobian_dxdr->resize(dim, dim);
        _prop->inv_jacobian_drdx->resize(dim, dim);
        (*_prop->shape_r) = .0;
        (*_prop->dshape_dr) = .0;
        (*_prop->dshape_dx) = .0;
        (*_prop->jacobian_dxdr) = .0;
        (*_prop->inv_jacobian_drdx) = .0;
        x.resize(dim, nnodes);
        MeshLib::IElementCoordinatesMapping *ele_local_coord = _ele->getMappedCoordinates();
        for (size_t i=0; i<nnodes; i++)
            for (size_t j=0; j<dim; j++)
                x(j,i) = ele_local_coord->getNodePoint(i)->getData()[j];
    }

    virtual const CoordMappingProperties* getProperties() const {return _prop;};

    /// compute mapping properties at the given location in natural coordinates
    virtual const CoordMappingProperties* compute(const double* natural_pt) 
    {
        //prepare
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();

        //shape, dshape
        MathLib::Matrix<double> *shape  = _prop->shape_r;
        MathLib::Matrix<double> *dshape_dr = _prop->dshape_dr;
        _shape->computeShapeFunction(natural_pt, (double*)shape->getData());
        _shape->computeGradShapeFunction(natural_pt, (double*)dshape_dr->getData());

        //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
        MathLib::Matrix<double> *jac = _prop->jacobian_dxdr;
        (*jac) = .0;
        for (size_t i_r=0; i_r<dim; i_r++) {
            for (size_t j_x=0; j_x<dim; j_x++) {
                for (size_t k=0; k<nnodes; k++) {
                    (*jac)(i_r,j_x) += (*dshape_dr)(i_r, k) * x(j_x,k);
                }
            }
        }

        //determinant of J
        double det_j = jac->determinant();
        _prop->det_jacobian = det_j;
        if (det_j>.0) {
            //inverse of J
            MathLib::Matrix<double> *inv_jac = _prop->inv_jacobian_drdx;
            jac->inverse(inv_jac);

            //dshape/dx = invJ * dNdr
            MathLib::Matrix<double> *dshape_dx = _prop->dshape_dx;
            (*dshape_dx) = .0;
            inv_jac->multiply(*dshape_dr, *dshape_dx);
        } else {
            std::cout << "***error: det_j is not positive." << std::endl;
        }

        return _prop;
    };

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param natural_pt
    /// @return physical_pt
    void mapToPhysicalCoordinates(const double* natural_pt, double* physical_pt)
    {
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();
        MathLib::Matrix<double> *shape  = _prop->shape_r;
        _shape->computeShapeFunction(natural_pt, (double*)shape->getData());

        for (size_t i=0; i<dim; i++) {
            physical_pt[i] = .0;
            for (size_t j=0; j<nnodes; j++)
                physical_pt[i] += (*shape)(0,j) + x(i,j);
        }
    }

    /// compute physical coordinates at the given natural coordinates
    /// \f[
    ///    \mathbf{x} = \mathbf{N(r)} * \mathbf{X}
    /// \f]
    ///
    /// @param natural_pt
    /// @return physical_pt
    void mapToPhysicalCoordinates(const CoordMappingProperties* prop, double* physical_pt)
    {
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();
        MathLib::Matrix<double> *shape  = prop->shape_r;

        for (size_t i=0; i<dim; i++) {
            physical_pt[i] = .0;
            for (size_t j=0; j<nnodes; j++)
                physical_pt[i] += (*shape)(0,j) * x(i,j);
        }
    }

    /// compute natural coordinates at the given natural coordinates.
    /// Assuming \f$ r=0 \f$ at \f$ x = \bar{x}_{avg} \f$, natural coordinates can be calculated as
    /// \f[
    ///    \mathbf{r} = (\mathbf{J}^-1)^T * (\mathbf{x} - \bar{\mathbf{x}}_{avg})
    /// \f]
    ///
    /// @param physical_pt
    /// @return natural_pt
    void mapFromPhysicalCoordinates(const double* physical_pt, double* natural_pt)
    {
        const size_t dim = _ele->getDimension();
        const size_t nnodes = _ele->getNumberOfNodes();
        MathLib::Matrix<double> *inv_jac = _prop->inv_jacobian_drdx;

        // calculate dx which is relative coordinates from element center
        std::vector<double> dx(dim, .0);
        // x_avg = sum_i {x_i} / n
        for (size_t i=0; i<dim; i++)
            for (size_t j=0; j<nnodes; j++)
                dx[i] += x(i,j);
        for (size_t i=0; i<dim; i++)
            dx[i] /= (double)nnodes;
        // dx = pt - x_avg
        for (size_t i=0; i<dim; i++)
            dx[i] = physical_pt[i] -dx[i];

        // r = invJ^T * dx
        for (size_t i=0; i<dim; i++) {
            natural_pt[i] = 0.0;
            for (size_t j=0; j<dim; j++)
                natural_pt[i] += (*inv_jac)(j * dim, i) * dx[j];
        }

    }
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

    virtual void computeMappingFunctions(double* natural_pt) 
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
