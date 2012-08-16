/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemSourceTerm.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/Tools.h"

#include "NumLib/Function/ITXFunction.h"

#include "IFemNeumannBC.h"

namespace SolutionLib
{

/**
 * 
 */
class FemSourceTerm : public IFemNeumannBC 
{
public:
    /// 
    FemSourceTerm(const MeshLib::IMesh *msh, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
    {
        _msh = msh;
        _geo = geo;
        _bc_func = func->clone();
        _is_transient = false;
        _do_setup = true;
    }

    /// setup  
    void setup();

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
    const MeshLib::IMesh* _msh;
    const GeoLib::GeoObject *_geo;
    NumLib::ITXFunction *_bc_func;
    std::vector<size_t> _vec_nodes;
    std::vector<double> _vec_values;
    bool _is_transient;
    bool _do_setup;

};

}
