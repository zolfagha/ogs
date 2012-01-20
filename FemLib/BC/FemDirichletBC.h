
#pragma once

#include "FemLib/FemFunction.h"

#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

namespace FemLib
{
/**
 * DirichletBC class
 */
template<typename Tval, typename Tpos, typename Tgeo>
class FemDirichletBC //per variable?
{
public:
    FemDirichletBC(FEMNodalFunction<Tval,Tpos> *fem, Tgeo *geo, MathLib::IFunction<Tval, Tpos> *func)
    {
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(fem->getMesh(), geo, &_vec_nodes);
        // set values
        _vec_values.reserve(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const Tpos* x = _vec_nodes[i]->getData();
            _vec_values[i] = func->eval(*x);
        }
    }

    template<typename Tmat, typename Tvec>
    void apply( Tmat* globalA, Tvec* globalRHS ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }


private:
    std::vector<MeshLib::Node<Tpos>*> _vec_nodes;
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
