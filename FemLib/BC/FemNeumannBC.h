
#pragma once

#include <map>
#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "NumLib/Function/Function.h"

#include "FemLib/Function/FemFunction.h"
#include "IFemBC.h"


namespace FemLib
{

class IFemNeumannBC
{
public:
    ///
    virtual ~IFemNeumannBC() {};

    /// setup B.C.
    virtual void setup() = 0;

    ///
    virtual std::vector<size_t>& getListOfBCNodes() = 0;
    virtual std::vector<double>& getListOfBCValues() = 0;
};

/**
 * Neumann BC
 */
template<typename Tval, typename Tflux>
class FemNeumannBC : public IFemNeumannBC //: IFemBC, public MathLib::TemplateSpatialFunction<Tflux>
{
public:
    /// 
    FemNeumannBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
    {
        _var = var;
        _geo = geo;
        _bc_func = func->clone();
        _is_transient = false;
        _do_setup = true;
    }

    /// setup BC. 
    void setup()
    {
        if (!_do_setup) return;
        if (!_is_transient) _do_setup = false;

        MeshLib::IMesh *msh = (MeshLib::IMesh*)_var->getMesh();
        LagrangianFeObjectContainer* feObjects = _var->getFeObjectContainer();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // distribute to RHS
        _vec_values.resize(_vec_nodes.size());
        if (msh->getDimension()==1) {
            // no need to integrate
            // get discrete values at nodes
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                const GeoLib::Point* x = msh->getNodeCoordinatesRef(_vec_nodes[i]);
                _bc_func->eval(*x, _vec_values[i]);
            }
        } else {
            // find edge elements on the geo
            std::vector<MeshLib::IElement*> vec_edge_eles;
            MeshLib::findBoundaryElementsOnGeometry(msh, _geo, &vec_edge_eles);
            // for each edge elements found
            std::map<size_t, Tval> map_nodeId2val;
            for (size_t i=0; i<vec_edge_eles.size(); i++) {
                MeshLib::IElement *e = vec_edge_eles[i];
                const size_t edge_nnodes = e->getNumberOfNodes();
                // set values at nodes
                LocalVector nodal_val(edge_nnodes);
                for (size_t i_nod=0; i_nod<edge_nnodes; i_nod++) {
                    const GeoLib::Point* x = msh->getNodeCoordinatesRef(e->getNodeID(i_nod));
                    double v = .0;
                    _bc_func->eval(*x, v);
                    nodal_val(i_nod) = v;
                } 
                // compute integrals
                IFiniteElement *fe_edge = feObjects->getFeObject(*e);
                fe_edge->getIntegrationMethod()->initialize(*e, 2);
                //IFiniteElement *fe_edge = _var->getFiniteElement(*e);
                LocalVector result(edge_nnodes);
                LocalMatrix M(edge_nnodes, edge_nnodes);
                M = .0;
                IFemNumericalIntegration *q = fe_edge->getIntegrationMethod();
                double x_ref[3];
                for (size_t j=0; j<q->getNumberOfSamplingPoints(); j++) {
                    q->getSamplingPoint(j, x_ref);
                    fe_edge->computeBasisFunctions(x_ref);
                	fe_edge->integrateWxN(j, 1., M);
                }
                result = M * nodal_val;
                //M.axpy(1.0, &nodal_val[0], 0.0, &result[0]);
                // add into RHS values
                for (size_t k=0; k<edge_nnodes; k++)
                    map_nodeId2val[e->getNodeID(k)] += result[k];
            }
            for (size_t i=0; i<_vec_nodes.size(); i++) {
                _vec_values[i] = map_nodeId2val[_vec_nodes[i]];
            }
        }
    }

//    ///
//    template<typename T>
//    void apply( T* globalRHS )
//    {
//        for (size_t i=0; i<this->getNumberOfConditions(); i++)
//            globalRHS[this->getConditionDoF(i)] -= this->getConditionValue(i);
//    }


    size_t getNumberOfConditions() const
    {
        return _vec_nodes.size();
    }

    size_t getConditionDoF( size_t i ) const
    {
        return _vec_nodes[i];
    }

    double getConditionValue( size_t i ) const
    {
        return _vec_values[i];
    }

//    void eval(const MathLib::SpatialPosition& x, Tflux &v)
//    {
//        _bc_func->eval(x, v);
//    }

    FemNeumannBC<Tval, Tflux>* clone() const
    {
        FemNeumannBC<Tval, Tflux> *f = new FemNeumannBC<Tval, Tflux>(_var, _geo, _is_transient, _bc_func);
        return f;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<Tval>& getListOfBCValues() {return _vec_values;};

private:
    // node id, var id, value
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<Tval> _vec_values;
    bool _is_transient;
    bool _do_setup;

};



}
