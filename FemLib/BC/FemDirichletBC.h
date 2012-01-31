
#pragma once

#include "MathLib/Function/Function.h"
#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/Function/FemFunction.h"
#include "IFemBC.h"

namespace FemLib
{
//----------------------------------------------------------

class IDirichletBCMethod
{
public:
//    virtual void apply(int linearEqs, DirichletBC &bc) = 0;
};

template<typename Tmat, typename Tvec, typename Tval>
class DiagonalizeMethod : public IDirichletBCMethod
{
private:
    void diagonalize(Tmat *eqsA, Tvec *eqsRHS, size_t id, Tval x);

public:
    void apply( Tmat* globalA, Tvec* globalRHS, std::vector<size_t> vec_nodes, std::vector<Tval> vec_values)
    {
        for (size_t i=0; i<vec_nodes.size(); i++)
            diagonalize(globalA, globalRHS, vec_nodes[i], vec_values[i]);
    }
};

template <> void DiagonalizeMethod<MathLib::Matrix<double>,std::vector<double>, double>::diagonalize(MathLib::Matrix<double> *eqsA, std::vector<double> *eqsRHS, size_t id, double x)
{
    const size_t n_cols = eqsA->getNCols();
    //A(k, j) = 0.
    for (size_t j=0; j<n_cols; j++)
        (*eqsA)(id, j) = .0;
    //b_i -= A(i,k)*val, i!=k
    for (size_t j=0; j<n_cols; j++)
        (*eqsRHS)[j] -= (*eqsA)(j, id)*x;
    //b_k = val
    (*eqsRHS)[id] = x;
    //A(i, k) = 0., i!=k
    for (size_t j=0; j<n_cols; j++)
      (*eqsA)(j, id) = .0;
    //A(k, k) = 1.0
    (*eqsA)(id, id) = 1.0; //=x
}


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


/**
 * DirichletBC class
 */
template<typename Tval>
class FemDirichletBC : IFemBC
{
public:
    ///
    FemDirichletBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, MathLib::IFunction<Tval, GeoLib::Point> *bc_func, IDirichletBCMethod *method)
    {
        _var = var;
        _geo = geo;
        _bc_func = bc_func;
        _method = method;
    }
    virtual ~FemDirichletBC()
    {
    }

    /// setup B.C.
    void setup()
    {
        const MeshLib::IMesh *msh = _var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // set values
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = msh->getNodeCoordinatesRef(_vec_nodes[i]);
            _vec_values[i] = _bc_func->eval(*x);
        }
    }

    /// apply B.C.
    template<typename Tmat, typename Tvec>
    void apply( Tmat* globalA, Tvec* globalRHS ) 
    {
        DiagonalizeMethod<Tmat, Tvec, Tval> method;
        method.apply(globalA, globalRHS, _vec_nodes, _vec_values);
    }


private:
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::IFunction<Tval, GeoLib::Point> *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<Tval> _vec_values;
    // node id, var id, value
    IDirichletBCMethod *_method;
};


}
