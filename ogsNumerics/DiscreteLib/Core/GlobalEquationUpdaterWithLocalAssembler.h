/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GlobalEquationUpdatorWithLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"
#include "IElemenetWiseLinearEquationLocalAssembler.h"

namespace DiscreteLib
{

template <class T_LOCAL>
class GlobalEquationUpdaterWithLocalAssembler
{
public:
    typedef T_LOCAL LocalAssemblerType;

    GlobalEquationUpdaterWithLocalAssembler(MeshLib::IMesh* msh, DofEquationIdTable* dofManager, LocalAssemblerType* a)
    : _msh(msh), _dofManager(dofManager), _e_assembler(a)
    {

    }

    void update(const MeshLib::IElement &e, MathLib::ILinearEquation &eqs)
    {
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap_row;
        std::vector<size_t> local_dofmap_column;
        LocalEquation localEQS;

        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        _dofManager->mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_column, local_dofmap_row); //TODO order
        //dofManager.mapEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
        // local assembly
        localEQS.create(local_dofmap_column.size());
        _e_assembler->assembly(e, localEQS);
        // update global
        eqs.addAsub(local_dofmap_row, local_dofmap_column, *localEQS.getA());
        eqs.addRHSsub(local_dofmap_row, localEQS.getRHS());
    }

private:
    MeshLib::IMesh* _msh;
    DofEquationIdTable* _dofManager;
    LocalAssemblerType* _e_assembler;
};

} //end
