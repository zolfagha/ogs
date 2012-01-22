
#pragma once

#include "FemLib/FemFunction.h"

#include "MathLib/Function/Function.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

namespace FemLib
{

/**
 * Neumann BC
 */
template<typename Tval, typename Tpos, typename Tgeo>
class FemNeumannBC
{
public:
    /// 
    FemNeumannBC(FEMNodalFunction<Tval,Tpos> *fem, Tgeo *geo, MathLib::IFunction<Tval, Tpos> *func) 
    {
        const MeshLib::IMesh<Tpos> *msh = fem->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, geo, &_vec_nodes);
        // get discrete values at nodes
        std::vector<Tval> vec_nod_values(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const Tpos* x = _vec_nodes[i]->getData();
            vec_nod_values[i] = func->eval(*x);
        }
        // distribute to nodes
        _vec_values.reserve(_vec_nodes.size());
        if (fem->getDimension()==1) {
          // no need to integrate
          _vec_values.assign(_vec_values.begin(), _vec_values.end());
        } else {
          // integrate over geo
          // find elements with boundaries on the geo
          std::vector<size_t> vec_ele;
          MeshLib::findElementsOnGeometry(msh, geo, &vec_ele);
          // for each elements found
          IFiniteElement *fe = fem->getFiniteElement();
          for (size_t i=0; i<vec_ele.size(); i++) {
            const MeshLib::IElement *e = msh->getElemenet(i);
            fe->configure(e);
            // find edges on the geo
            std::vector<size_t> vec_edges;
            // for each edges
            for (size_t j=0; j<vec_edges.size(); j++) {
              IFiniteElement *fe_bnd = fe->getBoundaryFE(j);
              // compute integrals
              std::vector<double> nodal_val(fe_bnd->getNumberOfDOFs());
              std::vector<double> result(fe_bnd->getNumberOfDOFs());
              fe_bnd->computeIntTestShapeNodalVal(&nodal_val[0], &result[0]);
              // add into nodal values
              const size_t edge_nodes = 0;
              std::vector<size_t> localNodeId2Global;
              for (size_t k=0; k<edge_nodes; k++)
                _vec_values[localNodeId2Global[k]] += result[k];
            }
          }
        }

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

    template<typename T>
    void apply( T* globalRHS ) 
    {
        for (size_t i=0; i<this->getNumberOfConditions(); i++)
            (*globalRHS)[this->getConditionDoF(i)] += this->getConditionValue(i);
    }



private:
    // node id, var id, value
    std::vector<MeshLib::Node<Tpos>*> _vec_nodes;
    std::vector<Tval> _vec_values;
};



}
