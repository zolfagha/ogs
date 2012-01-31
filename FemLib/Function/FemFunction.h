
#pragma once

#include "MathLib/Vector.h"
#include "MathLib/Function/Function.h"
#include "MeshLib/Core/IMesh.h"


#include "FemLib/Core/IFemElement.h"
#include "FemLib/FemElementObjectContainer.h"

namespace FemLib
{

struct PolynomialOrder
{
    enum type {
        Linear = 1,
        Quadratic = 2,
        INVALID = -1
    };
};

/**
 * \brief Template class for FEM node-based functions
 */
template<typename Tvalue>
class TemplateFEMNodalFunction : public MathLib::IFunction<Tvalue, GeoLib::Point>
{
public:
    TemplateFEMNodalFunction(MeshLib::IMesh* msh, PolynomialOrder::type order) {
        _msh = msh;
        _order = order;
        size_t nnodes = msh->getNumberOfNodes();
        _nodal_values = new Tvalue[nnodes];
    }

    virtual ~TemplateFEMNodalFunction() {
        delete _nodal_values;
    }

    size_t getDimension() const
    {
      return _msh->getCoordinateSystem()->getDimension();
    }

    const MeshLib::IMesh* getMesh() const {
        return _msh;
    }

    Tvalue eval(const GeoLib::Point &pt) {
        return _nodal_values[0];
    };

    Tvalue& getValue(int node_id) const {
        return _nodal_values[node_id];
    }

    IFiniteElement* getFiniteElement(MeshLib::IElement *e)
    {
        _feObjects.setPolynomialOrder(_order);
        IFiniteElement* fe = _feObjects.getFeObject(e);
        fe->configure(_msh, e);
        fe->getIntegrationMethod()->initialize(e, 2);

        return fe;
    }

    void setNodalValues( Tvalue* x ) 
    {
        for (size_t i=0; i<this->getNumberOfNodes(); i++)
            _nodal_values[i] = x[i];
    }

    size_t getNumberOfNodes() const 
    {
        return _msh->getNumberOfNodes();
    }

private:
    Tvalue* _nodal_values;
    //FEMInterpolation* _fe;
    MeshLib::IMesh* _msh;
    PolynomialOrder::type _order;
    LagrangianFeObjectContainer _feObjects;
};

typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar;
typedef TemplateFEMNodalFunction<MathLib::Vector2D> FemNodalFunctionVector2d;
typedef TemplateFEMNodalFunction<MathLib::Vector3D> FEMNodalFunctionVector3d;


/**
 * \brief Template class for FEM element-based functions
 */
template<typename Tvalue>
class TemplateFEMElementalFunction : public MathLib::IFunction<Tvalue,GeoLib::Point>
{
public:
    Tvalue eval(const GeoLib::Point &pt) {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};

/**
 * \brief Template class for FEM integration point-based functions
 */
template<typename Tvalue>
class TemplateFEMIntegrationPointFunction : public MathLib::IFunction<Tvalue,GeoLib::Point>
{
public:
    TemplateFEMIntegrationPointFunction(MeshLib::IMesh* msh) {
        _msh = msh;
        _values.resize(msh->getNumberOfElements());
    };

    Tvalue eval(const GeoLib::Point &pt) {
        throw std::exception("The method or operation is not implemented.");
    };

    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q ) 
    {
        assert(ip<_values[i_e].size());
        _values[i_e][ip] = q;
    }

    void setNumberOfIntegationPoints(size_t i_e, size_t n) {
        _values[i_e].resize(n);
    }

private:
    MeshLib::IMesh* _msh;
    std::vector<std::vector<Tvalue>> _values;
};

typedef TemplateFEMIntegrationPointFunction<double> FEMIntegrationPointFunctionScalar2d;
typedef TemplateFEMIntegrationPointFunction<MathLib::Vector2D> FEMIntegrationPointFunctionVector2d;

}



