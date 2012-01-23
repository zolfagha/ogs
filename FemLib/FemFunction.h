
#pragma once

#include "MathLib/Vector.h"
#include "MeshLib/Core/IMesh.h"
#include "NumLib/IFunction.h"
#include "FemElement.h"

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


template<typename Tvalue>
class TemplateFEMNodalFunction : public NumLib::IFunction<Tvalue>
{
public:
    TemplateFEMNodalFunction(MeshLib::IMesh* msh, PolynomialOrder::type order) {
        _msh = msh;
        _order = order;
    }

    size_t getDimension() const
    {
      return _msh->getDimension();
    }

    const MeshLib::IMesh* getMesh() const {
        return _msh;
    }

    Tvalue& getValue(GeoLib::Point &pt) const {
        return _nodal_values[0];
    };

    Tvalue& getValue(int node_id) const {
        return _nodal_values[node_id];
    }

    IFiniteElement* getFiniteElement() const
    {
        throw std::exception("The method or operation is not implemented.");
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
};

typedef TemplateFEMNodalFunction<double> FemNodalFunctionScalar2d;
typedef TemplateFEMNodalFunction<double> FEMNodalFunctionScalar3d;
typedef TemplateFEMNodalFunction<MathLib::Vector2D> FemNodalFunctionVector2d;
typedef TemplateFEMNodalFunction<MathLib::Vector3D> FEMNodalFunctionVector3d;

template<typename Tvalue>
class TemplateFEMElementalFunction : public NumLib::IFunction<Tvalue>
{
public:
    Tvalue& getValue(GeoLib::Point &pt) const {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};

template<typename Tvalue>
class TemplateFEMIntegrationPointFunction : public NumLib::IFunction<Tvalue>
{
public:
    TemplateFEMIntegrationPointFunction(MeshLib::IMesh* msh) {
        _msh = msh;
    };
    Tvalue& getValue(GeoLib::Point &pt) const {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
    void setIntegrationPointValue( size_t i_e, size_t ip, Tvalue &q ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

private:
    Tvalue* _ele_values;
    MeshLib::IMesh* _msh;
};

typedef TemplateFEMIntegrationPointFunction<double> FEMIntegrationPointFunctionScalar2d;
typedef TemplateFEMIntegrationPointFunction<MathLib::Vector2D> FEMIntegrationPointFunctionVector2d;

}



