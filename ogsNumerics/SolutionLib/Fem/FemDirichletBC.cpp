/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemDirichletBC.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemDirichletBC.h"

#include "FemLib/BC/DirichletBC2FEM.h"

namespace SolutionLib
{

/// setup B.C.
void FemDirichletBC::setup()
{
    if (!_do_setup) return;
    if (!_is_transient) _do_setup = false;

    _msh->setCurrentOrder(_order);
    FemLib::DirichletBC2FEM convert(*_msh, *_geo, *_bc_func, _vec_nodes, _vec_values);

    if (!_is_transient)
        _do_setup = false;
}

}
