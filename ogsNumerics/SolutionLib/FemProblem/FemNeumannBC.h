/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNeumannBC.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <map>
#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "NumLib/Function/Function.h"

#include "FemLib/Function/FemFunction.h"
#include "FemLib/BC/NeumannBC2FEM.h"

namespace SolutionLib
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
class FemNeumannBC : public IFemNeumannBC //: IFemBC, public MathLib::TemplateSpatialFunction<Tflux>
{
public:
    /// 
    FemNeumannBC(MeshLib::IMesh *msh, FemLib::LagrangianFeObjectContainer* feObjects, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
    {
        _msh = msh;
        _feObjects = feObjects;
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

        FemLib::NeumannBC2FEM convert(*_msh, *_feObjects, *_geo, *_bc_func, _vec_nodes, _vec_values);

        if (!_is_transient)
            _do_setup = false;
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

    FemNeumannBC* clone() const
    {
        FemNeumannBC *f = new FemNeumannBC(_msh, _feObjects, _geo, _bc_func);
        return f;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    // node id, var id, value
    MeshLib::IMesh* _msh;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    const GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;

};

/**
 * 
 */
class FemSourceTerm : public IFemNeumannBC 
{
public:
    /// 
    FemSourceTerm(MeshLib::IMesh *msh, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
    {
        _msh = msh;
        _geo = geo;
        _bc_func = func->clone();
        _is_transient = false;
        _do_setup = true;
    }

    /// setup  
    void setup()
    {
        if (!_do_setup) return;
        if (!_is_transient) _do_setup = false;

        MeshLib::findNodesOnGeometry(_msh, _geo, &_vec_nodes);
        _vec_values.resize(_vec_nodes.size());
        for (size_t i=0; i<_vec_nodes.size(); i++) {
            const GeoLib::Point* x = _msh->getNodeCoordinatesRef(_vec_nodes[i]);
            _bc_func->eval(x->getData(), _vec_values[i]);
        }

        if (!_is_transient)
            _do_setup = false;
    }

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


    FemSourceTerm* clone() const
    {
        FemSourceTerm *f = new FemSourceTerm(_msh, _geo, _bc_func);
        return f;
    }

    std::vector<size_t>& getListOfBCNodes() {return _vec_nodes;};
    std::vector<double>& getListOfBCValues() {return _vec_values;};

private:
    // node id, var id, value
    MeshLib::IMesh* _msh;
    const GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;

};

}
