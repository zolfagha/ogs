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

FemSourceTerm::FemSourceTerm(const MeshLib::IMesh *msh, const GeoLib::GeoObject *geo, NumLib::ITXFunction *func)
: _msh(msh), _geo(geo), _bc_func(func->clone())
{
    _is_transient = false;
    _do_setup = true;
}

FemSourceTerm::FemSourceTerm(const std::vector<size_t> &vec_node_id, const std::vector<double> &vec_node_values)
: _msh(NULL), _geo(NULL), _bc_func(NULL), _vec_nodes(vec_node_id), _vec_values(vec_node_values)
{
    _is_transient = false;
    _do_setup = false;
}

FemSourceTerm* FemSourceTerm::clone() const
{
    FemSourceTerm* f = NULL;
    if (_msh!=NULL)
        f = new FemSourceTerm(_msh, _geo, _bc_func);
    else
        f = new FemSourceTerm(_vec_nodes, _vec_values);

    return f;
}

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
