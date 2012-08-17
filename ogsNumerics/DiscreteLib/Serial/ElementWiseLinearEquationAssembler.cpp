/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseLinearEquationAssembler.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "ElementWiseLinearEquationAssembler.h"

#include <vector>

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/LocalDataType.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"


namespace DiscreteLib
{

void ElementWiseLinearEquationAssembler::assembly(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs)
{
    LocalEquation localEQS;
    std::vector<size_t> ele_node_ids, ele_node_size_order;
    std::vector<size_t> local_dofmap_row;
    std::vector<size_t> local_dofmap_column;
    const size_t n_ele = msh.getNumberOfElements();

    std::vector<double> local_u_n;
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        // get dof map
        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
        e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap_column, local_dofmap_row); //TODO order
        //dofManager.mapEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
        // local assembly
        localEQS.create(local_dofmap_column.size());
        _e_assembler->assembly(*e, localEQS);
        // update global
        eqs.addAsub(local_dofmap_row, local_dofmap_column, *localEQS.getA());
        eqs.addRHSsub(local_dofmap_row, localEQS.getRHS());
    }

    //apply ST
};

} //end

