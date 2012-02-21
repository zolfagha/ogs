
#pragma once

#include <vector>

#include "MathLib/Integration/GaussLegendre.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace FemLib
{

/**
 * \brief Interface for integration methods used in FEM
 */
class IFemNumericalIntegration
{
public:
    virtual void initialize(MeshLib::IElement &e, size_t n_sampl_level) = 0;
    virtual size_t getNumberOfSamplingPoints() const = 0;
    virtual void getSamplingPoint(size_t igp, double *) const = 0;
    virtual double getWeight(size_t igp) const = 0;
};

/**
 * \brief Dummy class for analytical integration methods which are not using sampling.
 */
class FemIntegrationAnalytical : public IFemNumericalIntegration
{
public:
    void initialize(MeshLib::IElement&, size_t) {};
    size_t getNumberOfSamplingPoints() const {return 1;};
    void getSamplingPoint(size_t, double*) const {};
    double getWeight(size_t igp) const {return 1.0;};
};

/**
 * \brief Gauss–Legendre integration class
 * 
 * 
 */
class FemIntegrationGaussBase : public IFemNumericalIntegration
{
public:
    FemIntegrationGaussBase()
    {
        _n_sampl_level = 0;
        _n_sampl_pt = 0;
    }

    void initialize(MeshLib::IElement &e, size_t n_sampl_level)
    {
        _n_sampl_level = n_sampl_level;
        _n_sampl_pt = getTotalNumberOfSamplingPoints(e, n_sampl_level);
    };

    size_t getSamplingLevel() const {return _n_sampl_level;};
    size_t getNumberOfSamplingPoints() const {return _n_sampl_pt;};

private:
    size_t _n_sampl_level;
    size_t _n_sampl_pt;

    virtual size_t getTotalNumberOfSamplingPoints(MeshLib::IElement &e, size_t n_sampl_level) const = 0;
    //size_t getTotalNumberOfSamplingPoints(MeshLib::IElement* e, size_t n_sampl_level) const
    //{
    //    switch (e->getElementType()) {
    //    case MeshLib::ElementType::LINE:
    //        return n_sampl_level;
    //    case MeshLib::ElementType::QUAD:
    //        return n_sampl_level*n_sampl_level;
    //    case MeshLib::ElementType::HEXAHEDRON:
    //        return n_sampl_level*n_sampl_level*n_sampl_level;
    //    case MeshLib::ElementType::TRIANGLE:
    //        if (n_sampl_level==1) return 1;
    //        else if (n_sampl_level==2) return 3;
    //        else return 6;
    //    }
    //    assert(false);
    //    return 0;
    //}
};

class FemIntegrationGaussLine : public FemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        x[0] = MathLib::GaussLegendre::getPoint(getSamplingLevel(), igp);
    }

    double getWeight(size_t igp) const
    {
        return MathLib::GaussLegendre::getWeight(getSamplingLevel(), igp);
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement &e, size_t n_sampl_level) const
    {
        return n_sampl_level;
    }
};

class FemIntegrationGaussQuad : public FemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        size_t nGauss = getSamplingLevel();
        size_t gp_r = igp/nGauss;
        size_t gp_s = igp%nGauss;
        x[0] = MathLib::GaussLegendre::getPoint(nGauss, gp_r);
        x[1] = MathLib::GaussLegendre::getPoint(nGauss, gp_s);
    }

    double getWeight(size_t igp) const
    {
        size_t nGauss = getSamplingLevel();
        size_t gp_r = igp/nGauss;
        size_t gp_s = igp%nGauss;
        return  MathLib::GaussLegendre::getWeight(nGauss, gp_r)*MathLib::GaussLegendre::getWeight(nGauss, gp_s);
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement &e, size_t n_sampl_level) const
    {
        return n_sampl_level*n_sampl_level;
    }
};

/**
 * \brief Gauss quadrature rule for triangles
 *
 * Gauss quadrature rule for triangles is originally given as
 * \f[
 *    \int F(x,y) dx dy = \int F(x(r, s), y(r, s)) j(r,s) dr ds \approx \frac{1}{2} \sum_i ( F(x(r_i, s_i), y(r_i, s_i)) w_i )
 * \f]
 *
 * To make it consistent with other elements, we rewrite the above formula as
 * \f[
 *    \int F(x,y) dx dy \approx \sum_i ( F(x(r_i, s_i), y(r_i, s_i)) w'_i )
 * \f]
 * by defining the new weight \f$ w'=\frac{1}{2} w \f$.
 */
class FemIntegrationGaussTriangle : public FemIntegrationGaussBase
{
public:
    void getSamplingPoint(size_t igp, double* x) const
    {
        switch (getNumberOfSamplingPoints())
        {
        case 1:
            break;
        case 3:
            getSamplePointTri3(igp, x);
            break;
        case 6:
            break;
        }
    }

    double getWeight(size_t igp) const
    {
        double w = .0;
        switch (getNumberOfSamplingPoints())
        {
        case 1:
            w= 1.0;
            break;
        case 3:
            w= 0.333333333333333; // = 1/3
        }

        return w*0.5;
    }

private:
    size_t getTotalNumberOfSamplingPoints(MeshLib::IElement &e, size_t n_sampl_level) const
    {
        if (n_sampl_level==1) return 1;
        else if (n_sampl_level==2) return 3;
        else return 6;
    }

    void getSamplePointTri1(size_t igp, double *pt) const
    {
        pt[0] = 0.333333333333333 ;
        pt[1] = 0.333333333333333 ;
    }

    void getSamplePointTri3(size_t igp, double *pt) const
    {
        switch(igp)
        {
        case 0:
            pt[0] = 0.166666666666667 ;
            pt[1] = 0.166666666666667 ;
            break;
        case 1:
            pt[0] = 0.666666666666667 ;
            pt[1] = 0.166666666666667 ;
            break;
        case 2:
            pt[0] = 0.166666666666667 ;
            pt[1] = 0.666666666666667 ;
            break;
        }
    }
};


class FemGaussIntegrationFactory
{
public:
    static FemIntegrationGaussBase* create(MeshLib::IElement *e) 
    {
        switch (e->getShapeType()) {
            case MeshLib::ElementShape::QUAD:
                return new FemIntegrationGaussQuad();
            case MeshLib::ElementShape::TRIANGLE:
                return new FemIntegrationGaussQuad();
            default:
                return 0;
        }
    }
};

}
