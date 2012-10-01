/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemNeumannBC.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemNeumannBC.h"

#include "logog.hpp"

#include "FemLib/BC/NeumannBC2FEM.h"

namespace SolutionLib
{

FemNeumannBC::FemNeumannBC(const MeshLib::IMesh *msh, FemLib::LagrangianFeObjectContainer* feObjects, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
: _msh(msh), _feObjects(feObjects), _geo(geo), _bc_func(func->clone())
{
    _is_transient = false;
    _do_setup = true;
}

FemNeumannBC::FemNeumannBC(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
: _msh(NULL), _feObjects(NULL), _geo(NULL), _bc_func(NULL), _vec_nodes(vec_node_id), _vec_values(vec_node_values)
{
    _is_transient = false;
    _do_setup = false;
}

FemNeumannBC* FemNeumannBC::clone() const
{
    FemNeumannBC* f = NULL;
    if (_msh!=NULL)
        f = new FemNeumannBC(_msh, _feObjects, _geo, _bc_func);
    else 
        f = new FemNeumannBC(_vec_nodes, _vec_values);

    return f;
}

/// setup BC.
void FemNeumannBC::setup(size_t order)
{
    if (!_do_setup) return;
    if (!_is_transient) _do_setup = false;

    _msh->setCurrentOrder(order);
    FemLib::NeumannBC2FEM convert(*_msh, *_feObjects, *_geo, *_bc_func, _vec_nodes, _vec_values);
    if (_vec_nodes.size()==0)
        INFO("***INFO: No Neumann BC found in FemDirichletBC::setup()");

    if (!_is_transient)
        _do_setup = false;
}

}
