
#pragma once

#include "MathLib/Function/Function.h"
#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/FemFunction.h"
#include "IFemBC.h"

namespace FemLib
{
/**
 * DirichletBC class
 */
template<typename Tval>
class FemDirichletBC : IFemBC
{
public:
    ///
    FemDirichletBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, MathLib::IFunction<Tval, GeoLib::Point> *bc_func)
    {
        _var = var;
        _geo = geo;
        _bc_func = bc_func;
    }

    /// setup B.C.
    void setup()
    {
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(_var->getMesh(), _geo, &_vec_nodes);
        // set values
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = _vec_nodes[i]->getData();
            _vec_values[i] = _bc_func->eval(*x);
        }
    }

    /// apply B.C.
    template<typename Tmat, typename Tvec>
    void apply( Tmat* globalA, Tvec* globalRHS ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }


private:
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::IFunction<Tval, GeoLib::Point> *_bc_func;
    std::vector<MeshLib::INode*> _vec_nodes;
    std::vector<Tval> _vec_values;
    // node id, var id, value
    //IDirichletBCMethod *_method;
};



//----------------------------------------------------------

class IDirichletBCMethod
{
public:
//    virtual void apply(int linearEqs, DirichletBC &bc) = 0;
};

class DiagonizeMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class ResizeMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class PenaltyMethod : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

class LagrangeMultiplier : public IDirichletBCMethod
{
public:
    void apply(int linearEqs, int &bc);
};

}
