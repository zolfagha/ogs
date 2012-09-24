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
