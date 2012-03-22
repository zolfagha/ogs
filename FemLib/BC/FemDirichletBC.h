
#pragma once

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
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

class DiagonalizeMethod : public IDirichletBCMethod
{
public:
    void apply( MathLib::ILinearEquations& eqs, const std::vector<size_t> &vec_nodes, const std::vector<double> &vec_values)
    {
        eqs.setKnownX(vec_nodes, vec_values);
    }
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


/**
 * DirichletBC class
 */
template<typename Tval>
class FemDirichletBC : IFemBC, public MathLib::IFunction<GeoLib::Point, Tval>
{
public:
    ///
    explicit FemDirichletBC(TemplateFEMNodalFunction<Tval> *var, GeoLib::GeoObject *geo, bool is_transient, MathLib::IFunction<GeoLib::Point, Tval> *bc_func, IDirichletBCMethod *method)
    {
        _var = var;
        _geo = geo;
        _bc_func = bc_func->clone();
        _method = method;
        _is_transient = is_transient;
        _do_setup = true;
    }

    virtual ~FemDirichletBC()
    {
    	//if (_method!=0) delete _method;
    	//_method = 0;
    }

    /// setup B.C.
    void setup()
    {
        if (!_do_setup) return;
        if (!_is_transient) _do_setup = false;

        const MeshLib::IMesh *msh = _var->getMesh();
        // pickup nodes on geo
        MeshLib::findNodesOnGeometry(msh, _geo, &_vec_nodes);
        // set values
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = msh->getNodeCoordinatesRef(_vec_nodes[i]);
            _bc_func->eval(*x, _vec_values[i]);
        }
        if (!_is_transient)
            _do_setup = false;
    }

    /// apply B.C.
    void apply( MathLib::ILinearEquations& eqs ) 
    {
        DiagonalizeMethod method;
        method.apply(eqs, _vec_nodes, _vec_values);
    }

    void eval(const GeoLib::Point &x, Tval &v)
    {
        _bc_func->eval(x, v);
    }

    MathLib::IFunction<GeoLib::Point,Tval>* clone() const
    {
        FemDirichletBC<Tval> *f = new FemDirichletBC<Tval>(_var, _geo, _is_transient, _bc_func, _method);
        return f;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<Tval>& getListOfBCValues() {return _vec_values;};

private:
    TemplateFEMNodalFunction<Tval> *_var;
    GeoLib::GeoObject *_geo;
    MathLib::IFunction<GeoLib::Point, Tval> *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<Tval> _vec_values;
    // node id, var id, value
    IDirichletBCMethod *_method;
    bool _is_transient;
    bool _do_setup;
};


}
