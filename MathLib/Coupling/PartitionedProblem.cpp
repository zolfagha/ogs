
#include "PartitionedProblem.h"

namespace MathLib
{

int PartitionedProblem::find(const ICoupledSystem& sub) const
{
    std::vector<ICoupledSystem*>::const_iterator itr = std::find(_list_subproblems.begin(), _list_subproblems.end(), &sub);
    if (itr!=_list_subproblems.end()) {
        return (itr - _list_subproblems.begin());
    } else {
        return -1;
    }
}

/// check consistency
bool PartitionedProblem::check() const
{
    bool flag = true;

    // check for subproblems
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        if (!_list_subproblems[i]->check())
            flag = false;
    }

    // check consistency in input and output of subproblems
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        const ICoupledSystem* subproblem = _list_subproblems[i];
        // check input parameters required for the subproblem
        const ParameterProblemMappingTable::ListOfInputVar &vec_registered_input_var = _map._list_subproblem_input_source[i];

        std::vector<size_t> list1(vec_registered_input_var.size());
        for (size_t j=0; j<list1.size(); j++) {
            list1[j] = vec_registered_input_var[j].first;
        }
        std::vector<size_t> list2(subproblem->getNumberOfInputParameters());
        for (size_t j=0; j<list2.size(); j++) {
            list2[j] = subproblem->getParameterIdForInput(j);
        }
        sort(list1.begin(), list1.end());
        sort(list2.begin(), list2.end());
        std::vector<size_t> list_diff(list1.size()+list2.size());
        std::vector<size_t>::iterator it = set_symmetric_difference(list1.begin(), list1.end(), list2.begin(), list2.end(), list_diff.begin());

        int n_diff = it-list_diff.begin();
        if (n_diff > 0) {
            for (int k=0; k<n_diff; k++) {
                std::cout << "*** Error: Inconsistent parameter found in subproblem " << i << " and param " << list_diff[k] << std::endl;
            }

            flag = false;
        }
    }

    return flag;
}

//input parameter
size_t PartitionedProblem::addParameter(const std::string &name)
{
    size_t n = _parameter_table.add(name);
    _map._map_paraId2subproblem.push_back(std::make_pair<ICoupledSystem*,size_t>(0, 0));
    _list_input_parameters.push_back(n);
    return n;
}

//output parameter
size_t PartitionedProblem::addParameter(const std::string &name, ICoupledSystem& sub_problem, size_t para_id_in_sub_problem)
{

    size_t var_id = _parameter_table.add(name);
    // set reference
    const Parameter* v = sub_problem.getOutput(para_id_in_sub_problem);
    _parameter_table.set(var_id, *v);
    // make a link between this and sub-problem variable
    ParameterProblemMappingTable::PairSysVarId parObj = std::make_pair(&sub_problem, para_id_in_sub_problem);
    _map._map_paraId2subproblem.push_back(parObj);
    // update a list of sub-problems
    if (find(sub_problem)<0) {
        addSubProblem(sub_problem);
    }

    return var_id;
}

void PartitionedProblem::connectInput(const std::string &this_para_name, ICoupledSystem &subproblem, size_t subproblem_para_id)
{
    size_t this_para_id = _parameter_table.find(this_para_name);
    int subproblem_id = find(subproblem);
    if (subproblem_id>=0) {
        _map._list_subproblem_input_source[subproblem_id].push_back(std::make_pair(subproblem_para_id, this_para_id));
    }
}

int PartitionedProblem::solve()
{
    return _algorithm->solve(_list_subproblems, _parameter_table, _map);
}

size_t PartitionedProblem::addSubProblem(ICoupledSystem &sub_problem)
{
    _list_subproblems.push_back(&sub_problem);
    _map._list_subproblem_input_source.resize(_list_subproblems.size());

    return _list_subproblems.size();
}

} //end
