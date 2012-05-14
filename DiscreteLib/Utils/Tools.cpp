
#include "Tools.h"


namespace DiscreteLib
{

/// create a subset of vector u corresponding to the given vector index
void getLocalVector2(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<IDiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u)
{
    local_u.clear();
    const size_t n_var = dofManager.getNumberOfVariables();
    for (size_t i=0; i<n_var; i++) {
        size_t var_order = 1; //TODO
        const size_t n_dof = list_vec_size_for_order[var_order-1];
        const DiscreteLib::IDiscreteVector<double> &var_u = *list_multiple_u[i];
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

/// create a subset of vector u corresponding to the given vector index
void getLocalVector(const DofEquationIdTable &dofManager, const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<IDiscreteVector<double>*> &list_multiple_u, LocalVector &local_u)
{
	std::vector<double> temp_v;
	getLocalVector2(dofManager, list_vec_entry_id, list_vec_size_for_order, list_multiple_u, temp_v);
	local_u.resize(temp_v.size());
	for (size_t i=0; i<temp_v.size(); i++)
		local_u[i] = temp_v[i];
}

} //end

