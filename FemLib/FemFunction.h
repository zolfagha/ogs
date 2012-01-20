
#pragma once

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


template<typename Tvalue, typename Tpos>
class FEMNodalFunction : public NumLib::IFunction<Tvalue,Tpos>
{
public:
    FEMNodalFunction(MeshLib::IMesh<Tpos>* msh, PolynomialOrder::type order) {
        _msh = msh;
        _order = order;
    }

    const MeshLib::IMesh<Tpos>* getMesh() const {
        return _msh;
    }

    Tvalue& getValue(Tpos &pt) const {
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
    MeshLib::IMesh<Tpos>* _msh;
    PolynomialOrder::type _order;
};


template<typename Tvalue, typename Tpos>
class FEMElementalFunction : public NumLib::IFunction<Tvalue,Tpos>
{
public:
    Tvalue& getValue(Tpos &pt) const {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};

template<typename Tvalue, typename Tpos>
class FEMIntegrationPointFunction : public NumLib::IFunction<Tvalue,Tpos>
{
public:
    FEMIntegrationPointFunction(MeshLib::IMesh<Tpos>* msh) {
        _msh = msh;
    };
    Tvalue& getValue(Tpos &pt) const {
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
    MeshLib::IMesh<Tpos>* _msh;
};


}

