
#include "AsyncPartitionedSystem.h"


namespace NumLib
{
double AsyncPartitionedSystem::suggestNext(const TimeStep &time_current) {
    double t;
    t = std::numeric_limits<double>::max();
    for (size_t i=0; i<_list_subproblems.size(); i++) {
    	MyTransientCoupledSystem *solution = _list_subproblems[i];
        double suggest = solution->suggestNext(time_current);
        if (suggest>.0)
            t = std::min(t, suggest);
    }
	return t;
}

void AsyncPartitionedSystem::getActiveProblems(const TimeStep &time, std::vector<MyCoupledSystem*> &list_active_problems)
{
    for (size_t i=0; i<_list_subproblems.size(); i++) {
    	MyTransientCoupledSystem *solution = _list_subproblems[i];
        if (solution->isAwake(time)) {
            list_active_problems.push_back(solution);
        }
    }
}

int AsyncPartitionedSystem::solveTimeStep(const TimeStep &time)
{

    std::cout << "->solve partitioned problems" << std::endl;

    // copy previous time step result to current one
    _vars_t_n.assign(_vars_t_n1);

    // list active problems
    std::vector<MyCoupledSystem*> list_active_problems;
    getActiveProblems(time, list_active_problems);

    if (list_active_problems.size()>0) {
        // solve
        for (size_t i=0; i<list_active_problems.size(); i++) {
        	MyTransientCoupledSystem *solution = (MyTransientCoupledSystem*)list_active_problems[i];
            solution->setCurrentTime(time);
        }

        //if (list_active_problems.size()==1) {
        //    list_active_problems[0]->solve();
        //} else {
            _algorithm->solve(list_active_problems, _vars_t_n, _vars_t_n1, _map);
        //}
    }


    return 0;
};

bool AsyncPartitionedSystem::isAwake(const TimeStep &time)
{
    for (size_t i=0; i<_list_subproblems.size(); i++) {
    	MyTransientCoupledSystem *solution = _list_subproblems[i];
        if (solution->isAwake(time))
            return true;
    }
    return false;
};

void AsyncPartitionedSystem::accept(const TimeStep &time)
{
    for (size_t i=0; i<_list_subproblems.size(); i++) {
    	MyTransientCoupledSystem *solution = _list_subproblems[i];
        solution->accept(time);
    }
}

int AsyncPartitionedSystem::find(const MyTransientCoupledSystem &sub) const
{
    typename std::vector<MyTransientCoupledSystem*>::const_iterator itr = std::find(_list_subproblems.begin(), _list_subproblems.end(), &sub);
    if (itr!=_list_subproblems.end()) {
        return (itr - _list_subproblems.begin());
    } else {
        return -1;
    }
}

/// check consistency
bool AsyncPartitionedSystem::check() const
{
    bool flag = true;

    // check for subproblems
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        if (!_list_subproblems[i]->check())
            flag = false;
    }

    // check consistency in input and output of subproblems
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        const MyTransientCoupledSystem* subproblem = _list_subproblems[i];
        // check input parameters required for the subproblem
        const typename MyVariableMappingTable::ListOfInputVar &vec_registered_input_var = _map._list_subproblem_input_source[i];

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

size_t AsyncPartitionedSystem::addParameter(const std::string &name)
{
    size_t n = _vars_t_n1.add(name);
    _map._map_paraId2subproblem.push_back(std::make_pair<MyTransientCoupledSystem*,size_t>(0, 0));
    _list_input_parameters.push_back(n);
    return n;
}

size_t AsyncPartitionedSystem::addParameter(const std::string &name, MyTransientCoupledSystem& sub_problem, size_t para_id_in_sub_problem)
{
    size_t var_id = _vars_t_n1.add(name);
    // set reference
    Variable* v = sub_problem.getParameter(para_id_in_sub_problem);
    _vars_t_n1.set(var_id, *v);
    // make a link between this and sub-problem variable
    _map._map_paraId2subproblem.push_back(std::make_pair(&sub_problem, para_id_in_sub_problem));
    // update a list of sub-problems
    if (find(sub_problem)<0) {
        addSubProblem(sub_problem);
    }

    return var_id;
}

void AsyncPartitionedSystem::connectInput(const std::string &this_para_name, MyTransientCoupledSystem &subproblem, size_t subproblem_para_id)
{
    int this_para_id = _vars_t_n1.find(this_para_name);
    if (this_para_id<0) return;

    int subproblem_id = find(subproblem);
    if (subproblem_id>=0) {
        _map._list_subproblem_input_source[subproblem_id].push_back(std::make_pair(subproblem_para_id, this_para_id));
    }
}

size_t AsyncPartitionedSystem::addSubProblem(MyTransientCoupledSystem &sub_problem)
{
    _list_subproblems.push_back(&sub_problem);
    _map._list_subproblem_input_source.resize(_list_subproblems.size());

    return _list_subproblems.size();
}
} //end
