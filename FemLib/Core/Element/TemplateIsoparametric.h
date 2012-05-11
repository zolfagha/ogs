#pragma once

#include <vector>
#include <cassert>

#include "MeshLib/Core/MeshGeometricProperties.h"

#include "FemLib/Core/Extrapolation/FemExtrapolation.h"
#include "FemLib/Core/CoordinatesMapping/FemNaturalCoordinates.h"
#include "TemplateFeBase.h"

namespace FemLib
{
             
/**
 * \brief Template class for any isoparametric FE classes
 *
 * \tparam T_FETYPE
 * \tparam N_VARIABLES
 * \tparam T_SHAPE
 * \tparam T_INTEGRAL
 * \tparam T_EXTRAPOLATE
 */
template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES, class T_SHAPE, class T_INTEGRAL, class T_EXTRAPOLATE>
class TemplateIsoparametric : public TemplateFeBase<T_FETYPE, N_VARIABLES>
{
public:
    ///
    TemplateIsoparametric(MeshLib::IMesh &msh) : TemplateFeBase<T_FETYPE, N_VARIABLES>(msh)
    {
        _mapping = new FemNaturalCoordinates(new T_SHAPE());
        _integration = new T_INTEGRAL();
        _is_basis_computed = false;
    };

    ///
    virtual ~TemplateIsoparametric() 
    {
        delete _mapping;
        delete _integration;
    };

    ///
    virtual IFemNumericalIntegration* getIntegrationMethod() const {return _integration;};

    /// initialize object for given mesh elements
    virtual void configure( MeshLib::IElement &e )
    {
    	TemplateFeBase<T_FETYPE, N_VARIABLES>::setElement(e);
        const MeshLib::IMesh* msh = TemplateFeBase<T_FETYPE, N_VARIABLES>::getMesh();
        if (e.getMappedCoordinates()==0) {
            MeshLib::IElementCoordinatesMapping *ele_map = 0;
            if (msh->getDimension() == e.getDimension()) {
                ele_map = new MeshLib::EleMapInvariant(*msh, e);
            } else {
                ele_map = new MeshLib::EleMapLocalCoordinates(*msh, e, *msh->getGeometricProperty()->getCoordinateSystem());
            }
            e.setMappedCoordinates(ele_map);
        }
        _mapping->initialize(e);
        _integration->initialize(e, 2);
        _is_basis_computed = false;
    }

    virtual void computeBasisFunctions(const double *x)
    {
        _mapping->compute(x);
        _is_basis_computed = true;
    }

    /// compute real coordinates from the given position in reference coordinates
    virtual void getRealCoordinates(double* x_real)
    {
    	assert(_is_basis_computed);
    	_mapping->mapToPhysicalCoordinates(_mapping->getProperties(), x_real);
    }

    virtual LocalMatrix* getBasisFunction()
    {
    	assert(_is_basis_computed);
        return _mapping->getProperties()->shape_r;
    }

    virtual LocalMatrix* getGradBasisFunction()
    {
    	assert(_is_basis_computed);
        return _mapping->getProperties()->dshape_dx;
    }

    virtual double getDetJ() const
    {
    	assert(_is_basis_computed);
        return _mapping->getProperties()->det_jacobian;
    }

    /// make interpolation from nodal values
    virtual double interpolate(double *natural_pt, double *nodal_values)
    {
        const CoordinateMappingProperty *prop = _mapping->compute(natural_pt);
        double *N = &(*prop->shape_r)(0,0);
        double v = .0;
        for (size_t i=0; i<TemplateFeBase<T_FETYPE, N_VARIABLES>::getNumberOfVariables(); i++)
            v+=N[i]*nodal_values[i];
        return v;
    }

//    /// compute an matrix M = Int{W^T F N} dV
//    virtual void integrateWxN(MathLib::SpatialFunctionScalar* f, LocalMatrix &mat)
//    {
//        const size_t n_gp = _integration->getNumberOfSamplingPoints();
//        double x[3], ox[3];
//        for (size_t i=0; i<n_gp; i++) {
//            _integration->getSamplingPoint(i, x);
//            _mapping->compute(x);
//            _mapping->mapToPhysicalCoordinates(_mapping->getProperties(), ox);
//            double fac = 1.0;
//            if (f!=0) {
//            	double v;
//                f->eval(ox, v);
//                fac *= v;
//            }
//            integrateWxN(i, fac, mat);
//        }
//    }
//
//    /// compute an matrix M = Int{W^T F dN} dV
//    virtual void integrateWxDN(MathLib::SpatialFunctionVector* f, LocalMatrix &mat)
//    {
//        const size_t n_gp = _integration->getNumberOfSamplingPoints();
//        double x[3], ox[3];
//        for (size_t i=0; i<n_gp; i++) {
//            _integration->getSamplingPoint(i, x);
//            _mapping->compute(x);
//            _mapping->mapToPhysicalCoordinates(_mapping->getProperties(), ox);
//            MathLib::Vector fac;
//            if (f!=0)
//                f->eval(ox, fac);
//            integrateWxDN(i, fac, mat);
//        }
//    }
//
//    /// compute an matrix M = Int{dW^T F dN} dV
//    virtual void integrateDWxDN(MathLib::SpatialFunctionScalar *f, LocalMatrix &mat)
//    {
//        const size_t n_gp = _integration->getNumberOfSamplingPoints();
//        double x[3], ox[3];
//        for (size_t i=0; i<n_gp; i++) {
//            _integration->getSamplingPoint(i, x);
//            const CoordinateMappingProperty *coord_prop = _mapping->compute(x);
//            _mapping->mapToPhysicalCoordinates(coord_prop,  ox);
//            double fac = 1.0;
//            if (f!=0) {
//            	double v;
//                f->eval(ox, v);
//                fac *= v;
//            }
//            integrateDWxDN(i, fac, mat);
//        }
//    }

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(size_t igp, double f, LocalMatrix &mat)
    {
    	assert(_is_basis_computed);
        const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
        LocalMatrix *basis = coord_prop->shape_r;
        LocalMatrix *test = coord_prop->shape_r;
        double fac = f * coord_prop->det_jacobian * _integration->getWeight(igp);
        //test->transposeAndMultiply(*basis, mat, fac);
        mat.noalias() += test->transpose() * (*basis) * fac;
    }

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(size_t igp, LocalMatrix &f, LocalMatrix &mat)
    {
    	assert(_is_basis_computed);
        const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
        LocalMatrix *dbasis = coord_prop->dshape_dx;
        LocalMatrix *test = coord_prop->shape_r;
        double fac = coord_prop->det_jacobian * _integration->getWeight(igp);
//        test->transposeAndMultiply(*dbasis, &f[0], mat, fac);
        mat.noalias() += test->transpose() * f * (*dbasis) * fac;
    }

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(size_t igp, LocalMatrix &f, LocalMatrix &mat)
    {
    	assert(_is_basis_computed);
        const CoordinateMappingProperty *coord_prop = _mapping->getProperties();
        LocalMatrix *dbasis = coord_prop->dshape_dx;
        LocalMatrix *dtest = coord_prop->dshape_dx;
        double fac = coord_prop->det_jacobian * _integration->getWeight(igp);
//        dtest->transposeAndMultiply(*dbasis, mat, fac);
        mat.noalias() += dtest->transpose() * f * (*dbasis) * fac;
    }

    void extrapolate(const std::vector<LocalVector> &gp_values, std::vector<LocalVector> &nodal_values)
    {
        T_EXTRAPOLATE extrapo;
        extrapo.extrapolate(*this, gp_values, nodal_values);
    }

private:
    FemNaturalCoordinates *_mapping;
    IFemNumericalIntegration* _integration;
    bool _is_basis_computed;
};

} //end
