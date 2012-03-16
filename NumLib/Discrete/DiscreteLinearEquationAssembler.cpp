
#pragma once

#include "DiscreteLinearEquationAssembler.h"

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"


namespace NumLib
{

/// create a subset of vector u corresponding to the given vector index
void getLocalVector(DiscreteLib::DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<DiscreteLib::DiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u)
{
    local_u.clear();
    const size_t n_var = dofManager.getNumberOfVariables();
    for (size_t i=0; i<n_var; i++) {
        size_t var_order = 1; //TODO
        const size_t n_dof = list_vec_size_for_order[var_order-1];
        const DiscreteLib::DiscreteVector<double> &var_u = *list_multiple_u[i];
        for (size_t j=0; j<n_dof; j++) {
            local_u.push_back(var_u[list_vec_entry_id[j]]);
        }
    }
    //const size_t n_dof = getNumberOfVariables();
    //for (size_t i=0; i<n_dof; i++) {
    //    const DofMap *dofMap = getVariableDoF(i);
    //    const size_t vec_size = list_vec_size_for_order[dofMap->getOrder()-1];
    //    const DiscreteVector<double> &dof_u = *list_multiple_u[i];
    //    for (size_t j=0; j<vec_size; j++) {
    //        local_u.push_back(dof_u[list_vec_entry_id[j]]);
    //    }
    //}
}


void ElementBasedTransientAssembler::assembly(MeshLib::IMesh &msh, DiscreteLib::DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs)
{
    const TimeStep &time = *_timestep;
    MathLib::DenseLinearEquations localEQS;
    std::vector<size_t> ele_node_ids, ele_node_size_order;
    std::vector<long> local_dofmap;
    const size_t n_ele = msh.getNumberOfElements();

    std::vector<double> local_u_n;
    for (size_t i=0; i<n_ele; i++) {
        MeshLib::IElement *e = msh.getElemenet(i);
        // get dof map
        e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);
        e->getListOfNumberOfNodesForAllOrders(ele_node_size_order);
        dofManager.mapEqsID(msh.getID(), ele_node_ids, local_dofmap);
        //dofManager.mapEqsID(ele_node_ids, ele_node_size_order, local_dofmap);
        // get previous time step results
        getLocalVector(dofManager, ele_node_ids, ele_node_size_order, *_u0, local_u_n);
        // local assembly
        localEQS.create(local_dofmap.size());
        _transient_e_assembler->assembly(time, *e, local_u_n, localEQS);
        // update global
        eqs.addAsub(local_dofmap, *localEQS.getA());
        eqs.addRHSsub(local_dofmap, localEQS.getRHS());
    }

    //apply ST
}


} //end

