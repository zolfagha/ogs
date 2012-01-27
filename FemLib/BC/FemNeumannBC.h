
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "FemLib/FemFunction.h"
#include "IFemBC.h"

namespace FemLib
{

/**
 * Neumann BC
 */
template<typename Tval>
class FemNeumannBC : IFemBC
{
public:
    /// 
    FemNeumannBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, MathLib::IFunction<Tval, GeoLib::Point> *func) 
    {
        _var = var;
        _geo = geo;
        _bc_func = func;
    }

    /// setup BC. 
    void setup()
    {
        MeshLib::IMesh *msh = (MeshLib::IMesh*)_var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // get discrete values at nodes
        std::vector<Tval> vec_nod_values(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = _vec_nodes[i]->getData();
            vec_nod_values[i] = _bc_func->eval(*x);
        }
        // distribute to nodes
        _vec_values.resize(_vec_nodes.size());
        if (_var->getDimension()==1) {
            // no need to integrate
            _vec_values.assign(_vec_values.begin(), _vec_values.end());
        } else {
            // integrate over geo
            // find edge elements on the geo
            std::vector<MeshLib::IElement*> vec_ele;
            MeshLib::findBoundaryElementsOnGeometry(msh, _geo, &vec_ele);
            // for each edge elements found
            for (size_t i=0; i<vec_ele.size(); i++) {
                MeshLib::IElement *e = vec_ele[i];
                IFiniteElement *fe_edge = _var->getFiniteElement(e);
                // compute integrals
                std::vector<double> nodal_val(fe_edge->getNumberOfVariables());
                std::vector<double> result(fe_edge->getNumberOfVariables());
                MathLib::Matrix<double> M(e->getNumberOfNodes(), e->getNumberOfNodes());
                fe_edge->integrateWxN(0, M);
                M.axpy(1.0, &nodal_val[0], 0.0, &result[0]);
                // add into nodal values
                const size_t edge_nodes = 0;
                std::vector<size_t> localNodeId2Global;
                for (size_t k=0; k<edge_nodes; k++)
                    _vec_values[localNodeId2Global[k]] += result[k];
            }
        }
    }

    ///
    template<typename T>
    void apply( T* globalRHS ) 
    {
        for (size_t i=0; i<this->getNumberOfConditions(); i++)
            (*globalRHS)[this->getConditionDoF(i)] += this->getConditionValue(i);
    }


    size_t getNumberOfConditions() const
    {
        throw std::exception("The method or operation is not implemented.");
    }

    size_t getConditionDoF( size_t i ) const
    {
        throw std::exception("The method or operation is not implemented.");
    }

    double getConditionValue( size_t i ) const
    {
        throw std::exception("The method or operation is not implemented.");
    }

private:
    // node id, var id, value
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::IFunction<Tval, GeoLib::Point> *_bc_func;
    std::vector<MeshLib::INode*> _vec_nodes;
    std::vector<Tval> _vec_values;
};



}
