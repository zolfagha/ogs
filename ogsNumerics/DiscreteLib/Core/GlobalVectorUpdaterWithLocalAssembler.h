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

template <class T_VALUE, class T_LOCAL>
class GlobalVectorUpdaterWithLocalAssembler
{
public:
    typedef T_LOCAL LocalAssemblerType;
    typedef IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;

    GlobalVectorUpdaterWithLocalAssembler(MeshLib::IMesh* msh, DofEquationIdTable* dofManager, LocalAssemblerType* a)
    : _msh(msh), _dofManager(dofManager), _e_assembler(a)
    {

    }

    void update(const MeshLib::IElement &e, GlobalVectorType &globalVec)
    {
        LocalVector localVec;
        std::vector<size_t> ele_node_ids, ele_node_size_order;
        std::vector<long> local_dofmap_row, local_dofmap_column;

        std::vector<T_VALUE> local_u_n;
        // get dof map
        e.getNodeIDList(e.getMaximumOrder(), ele_node_ids);
        e.getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        _dofManager->mapEqsID(_msh->getID(), ele_node_ids, local_dofmap_column, local_dofmap_row); //TODO order
        // local assembly
        localVec.resize(local_dofmap_column.size(), .0);
        _e_assembler->assembly(*e, localVec);
        // update global
        globalVec.addSubvector(local_dofmap_row, &localVec[0]);
    }

private:
    MeshLib::IMesh* _msh;
    DofEquationIdTable* _dofManager;
    LocalAssemblerType* _e_assembler;
};

} //end
