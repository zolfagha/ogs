
#pragma once

#include "MeshLib/Core/IMesh.h"
#include "NumLib/IFunction.h"
#include "Projection.h"

namespace FemLib
{

void mapNodeFunction2Element();
void mapElementFunction2Node();

struct LagrangeOrder
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
    FEMNodalFunction(MeshLib::IMesh* msh, LagrangeOrder::type order) {
        _msh = msh;
        _order = order;
    }

    Tvalue& getValue(Tpos &pt) {
        return _nodal_values[0];
    };

    Tvalue& getValue(int node_id) {
        return _nodal_values[node_id];
    }

    IFiniteElement* getFiniteElement() 
    {
        throw std::exception("The method or operation is not implemented.");
    }

private:
    Tvalue* _nodal_values;
    FEMInterpolation* _fe;
    MeshLib::IMesh* _msh;
    LagrangeOrder::type _order;
};


template<typename Tvalue, typename Tpos>
class FEMElementalFunction : public NumLib::IFunction<Tvalue,Tpos>
{
public:
    Tvalue& getValue(Tpos &pt) {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};

template<typename Tvalue, typename Tpos>
class FEMIntegrationPointFunction : public NumLib::IFunction<Tvalue,Tpos>
{
public:
    FEMIntegrationPointFunction(MeshLib::IMesh* msh) {
        _msh = msh;
    };
    Tvalue& getValue(Tpos &pt) {
        return _ele_values[0];
    };
    Tvalue& getValue(size_t id) {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
    MeshLib::IMesh* _msh;
};



}

