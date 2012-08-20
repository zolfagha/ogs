/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseLinearEquationUpdater.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "MeshLib/Core/IElement.h"
#include "DiscreteLib/Utils/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-wise global equation updater
 *
 * \tparam T_LOCAL Local assembler
 */
template <class T_LOCAL>
class ElementWiseLinearEquationUpdater
{
public:
    typedef T_LOCAL LocalAssemblerType;

    /**
     *
     * @param msh
     * @param dofManager
     * @param local_assembler
     */
    ElementWiseLinearEquationUpdater(MeshLib::IMesh* msh, DofEquationIdTable* dofManager, LocalAssemblerType* local_assembler)
    : _msh(msh), _dofManager(dofManager), _e_assembler(local_assembler)
    {

    }

    /**
     *
     * @param e     mesh element
     * @param eqs   global equation
     */
    void update(const MeshLib::IElement &e, MathLib::ILinearEquation &eqs)
    {
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<size_t> local_dofmap_row;
        //std::vector<size_t> local_dofmap_column;
        LocalEquation localEQS;

        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        _dofManager->mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_row); //TODO order
        //dofManager.mapEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
        // local assembly
        localEQS.create(local_dofmap_row.size());
        _e_assembler->assembly(e, localEQS);
        // update global
        eqs.addAsub(local_dofmap_row, local_dofmap_row, *localEQS.getA());
        eqs.addRHSsub(local_dofmap_row, localEQS.getRHS());
    }

private:
    MeshLib::IMesh* _msh;
    DofEquationIdTable* _dofManager;
    LocalAssemblerType* _e_assembler;
};

} //end
