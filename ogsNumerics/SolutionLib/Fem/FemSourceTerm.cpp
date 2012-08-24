/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemSourceTerm.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemSourceTerm.h"

#include "MeshLib/Tools/Tools.h"


namespace SolutionLib
{

/// setup
void FemSourceTerm::setup(size_t order)
{
    if (!_do_setup) return;
    if (!_is_transient) _do_setup = false;

    _msh->setCurrentOrder(order);
    MeshLib::findNodesOnGeometry(_msh, _geo, &_vec_nodes);
    _vec_values.resize(_vec_nodes.size());
    for (size_t i=0; i<_vec_nodes.size(); i++) {
        const GeoLib::Point* x = _msh->getNodeCoordinatesRef(_vec_nodes[i]);
        _bc_func->eval(x->getData(), _vec_values[i]);
    }

    if (!_is_transient)
        _do_setup = false;
}

}
