
#pragma once

#include <vector>
#include "DiscreteLib/Core/DiscreteVector.h"
#include "DiscreteLib/EquationId/DofEquationIdTable.h"

namespace DiscreteLib
{

/// create a subset of vector u corresponding to the given vector index
void getLocalVector(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<DiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u)
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

} //end

