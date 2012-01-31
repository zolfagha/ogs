#pragma once

#include "FemLib/Core/IFemElement.h"

namespace FemLib
{

/**
 * \brief Base for any isoparametric FE classes
 */
template <FiniteElementType::type T_FETYPE, size_t N_VARIABLES, class T_SHAPE, class T_INTEGRAL>
class FeBaseIsoparametric : public TemplateFeBase<T_FETYPE, N_VARIABLES>
{
public:
    FeBaseIsoparametric()
    {
        _mapping = new FemNaturalCoordinates(new T_SHAPE());
        _integration = new T_INTEGRAL();
    };

    virtual ~FeBaseIsoparametric() 
    {
        delete _mapping;
        delete _integration;
    };

    virtual IFemNumericalIntegration* getIntegrationMethod() const {return _integration;};

    /// initialize object for given mesh elements
    virtual void configure( MeshLib::IMesh * msh, MeshLib::IElement * e )
    {
        setElement(e);
        if (e->getMappedCoordinates()==0) {
            MeshLib::IElementCoordinatesMapping *ele_map = 0;
            if (msh->getCoordinateSystem()->getDimension() == e->getDimension()) {
                ele_map = new MeshLib::EleMapInvariant(msh, e);
            } else {
                ele_map = new MeshLib::EleMapLocalCoordinates(msh, e, msh->getCoordinateSystem());
            }
            e->setMappedCoordinates(ele_map);
        }
        _mapping->initialize(e);
    }

    virtual void computeBasisFunctions(const double *x)
    {
        _mapping->compute(x);
    }

    virtual MathLib::Matrix<double>* getBasisFunction()
    {
        return _mapping->getProperties()->shape_r;
    }

    virtual MathLib::Matrix<double>* getGradBasisFunction()
    {
        return _mapping->getProperties()->dshape_dx;
    }


    /// make interpolation from nodal values
    virtual double interpolate(double *natural_pt, double *nodal_values)
    {
        const CoordMappingProperties *prop = _mapping->compute(natural_pt);
        double *N = (double*)prop->shape_r->getData();
        double v = .0;
        for (size_t i=0; i<getNumberOfVariables(); i++)
            v+=N[i]*nodal_values[i];
        return v;
    }

    /// compute an matrix M = Int{W^T F N} dV
    virtual void integrateWxN(MathLib::IFunction<double, double*>* f, MathLib::Matrix<double> &mat)
    {
        const size_t n_gp = _integration->getNumberOfSamplingPoints();
        double x[3];
        for (size_t i=0; i<n_gp; i++) {
            _integration->getSamplingPoint(i, x);
            const CoordMappingProperties *coord_prop = _mapping->compute(x);
            MathLib::Matrix<double> *basis = coord_prop->shape_r;
            MathLib::Matrix<double> *test = coord_prop->shape_r;
            double fac = coord_prop->det_jacobian * _integration->getWeight(i);
            if (f!=0)
                fac *= f->eval(x);
            test->transposeAndMultiply(*basis, mat, fac);
        }
    }

    /// compute an matrix M = Int{W^T F dN} dV
    virtual void integrateWxDN(MathLib::IFunction<double*, double*>* f, MathLib::Matrix<double> &mat)
    {
        const size_t n_gp = _integration->getNumberOfSamplingPoints();
        double x[3];
        for (size_t i=0; i<n_gp; i++) {
            _integration->getSamplingPoint(i, x);
            const CoordMappingProperties *coord_prop = _mapping->compute(x);
            MathLib::Matrix<double> *dbasis = coord_prop->dshape_dx;
            MathLib::Matrix<double> *test = coord_prop->dshape_dx;
            double fac = coord_prop->det_jacobian * _integration->getWeight(i);
            if (f==0) {
                test->transposeAndMultiply(*dbasis, mat, fac);
            } else {
                double *v = f->eval(x);
                test->transposeAndMultiply(*dbasis, v, mat, fac);
            }
        }
    }

    /// compute an matrix M = Int{dW^T F dN} dV
    virtual void integrateDWxDN(MathLib::IFunction<double, double*> *f, MathLib::Matrix<double> &mat)
    {
        const size_t n_gp = _integration->getNumberOfSamplingPoints();
        double x[3];
        for (size_t i=0; i<n_gp; i++) {
            _integration->getSamplingPoint(i, x);
            const CoordMappingProperties *coord_prop = _mapping->compute(x);
            MathLib::Matrix<double> *dbasis = coord_prop->dshape_dx;
            MathLib::Matrix<double> *dtest = coord_prop->dshape_dx;
            double fac = coord_prop->det_jacobian * _integration->getWeight(i);
            if (f!=0)
                fac *= f->eval(x);
            dtest->transposeAndMultiply(*dbasis, mat, fac);
        }
    }

private:
    FemNaturalCoordinates *_mapping;
    IFemNumericalIntegration* _integration;

//protected:
//    virtual IFemShapeFunction* createShapeFunction() const = 0;
//    virtual IFemIntegration* createIntegrationMethod() const = 0;
};


typedef FeBaseIsoparametric<FiniteElementType::LINE2, 2, FemShapeLine2, FemIntegrationGaussLine> LINE2;
typedef FeBaseIsoparametric<FiniteElementType::LINE3, 3, FemShapeLine3, FemIntegrationGaussLine> LINE3;
typedef FeBaseIsoparametric<FiniteElementType::QUAD4, 4, FemShapeQuad4, FemIntegrationGaussQuad> QUAD4;
typedef FeBaseIsoparametric<FiniteElementType::QUAD8, 8, FemShapeQuad8, FemIntegrationGaussQuad> QUAD8;
typedef FeBaseIsoparametric<FiniteElementType::QUAD9, 9, FemShapeQuad9, FemIntegrationGaussQuad> QUAD9;
typedef FeBaseIsoparametric<FiniteElementType::TRI3, 3, FemShapeTriangle3, FemIntegrationGaussTriangle> TRI3;
typedef FeBaseIsoparametric<FiniteElementType::TRI6, 6, FemShapeTriangle6, FemIntegrationGaussTriangle> TRI6;

}
