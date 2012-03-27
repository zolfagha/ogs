
#pragma once

#include <cstddef>
#include <iostream>
#include <cmath>

#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/VariableContainer.h"
#include "NumLib/Coupling/VariableMappingTable.h"
#include "PartitionedAlgorithm.h"

namespace NumLib
{

/**
 * \brief Abstract class for iterative partitioned methods
 */
template <class T_CONVERGENCE_CHECK>
class AbstractIterativePartitionedMethod : public IPartitionedAlgorithm
{
	typedef MathLib::IFunction Variable;
	typedef VariableContainer MyNamedVariableContainer;
	typedef ICoupledSystem MyCoupledSystem;
	typedef VariableMappingTable MyVariableMappingTable;
public:
    AbstractIterativePartitionedMethod() : _max_itr(100), _epsilon(1e-3) {}

    AbstractIterativePartitionedMethod(double epsilon, size_t max_count) : _max_itr(max_count), _epsilon(epsilon) {}

    virtual ~AbstractIterativePartitionedMethod() {};

    /// get max. iteration
    size_t getMaximumIterationCounts() const {return _max_itr;};

    /// get epsilon
    double getEpsilon() const {return _epsilon;};

    /// get iteration count
    size_t getIterationCounts() const {return _itr_count;};

    /// solve
    int solve(std::vector<MyCoupledSystem*> &subproblems, MyNamedVariableContainer &vars, MyVariableMappingTable &mapping);

protected:
    //// check if solution is converged
    //bool isConverged(MyNamedVariableContainer& vars_prev, MyNamedVariableContainer& vars_current, double &v_diff);

    virtual void doPostAfterSolve( MyCoupledSystem& /*solution*/, MyNamedVariableContainer& /*vars*/, MyVariableMappingTable& /*mapping*/ )  {}

    virtual void doPostAfterSolveAll( MyNamedVariableContainer& /*vars*/, MyVariableMappingTable& /*mapping*/ ) {}

private:
    size_t _max_itr;
    size_t _itr_count;
    double _epsilon;
};

template <class T_CONVERGENCE_CHECK>
int AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>::solve(std::vector<MyCoupledSystem*> &subproblems, MyNamedVariableContainer &vars, MyVariableMappingTable &mapping)
{
    const size_t n_subproblems = subproblems.size();

    // initialize variables
    const size_t n_para = vars.size();
    for (size_t i=0; i<n_para; i++) {
        typename MyVariableMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
        MyCoupledSystem *solution = sysvarid.first;
        size_t internalId = sysvarid.second;
        if (solution!=0)
            vars.set(i, *const_cast<Variable*>(solution->getParameter(internalId)));
    }

    // obj to store previous state
    MyNamedVariableContainer vars_prev;

    // iteration start
    const size_t max_itr = getMaximumIterationCounts();
    bool is_converged = false;
    size_t i_itr = 0;
    double v_diff = .0;
    do {
        // copy current to prev
        vars.clone(vars_prev);
        // compute each solution
        for (size_t i=0; i<n_subproblems; i++) {
        	MyCoupledSystem *solution = subproblems[i];
            std::vector<MyVariableMappingTable::PairInputVar> &list_input_var = mapping._list_subproblem_input_source[i];
            // set input
            for (size_t j=0; j<list_input_var.size(); j++) {
                size_t local_var_id = list_input_var[j].first;
                size_t shared_var_id = list_input_var[j].second;
                solution->setParameter(local_var_id, vars.get(shared_var_id));
            }
            solution->solve();

            //
            doPostAfterSolve(*solution, vars, mapping);
        }
        doPostAfterSolveAll(vars, mapping);
        if (n_subproblems>1) {
            // check convergence
        	T_CONVERGENCE_CHECK check;
            is_converged = check.isConverged(vars_prev, vars, getEpsilon(), v_diff);
            std::cout << v_diff << " ";
        } else {
            is_converged = true;
        }
        ++i_itr;
    } while (!is_converged && i_itr<max_itr);
    std::cout << std::endl;
    std::cout << "iteration count=" << i_itr << ", error=" << v_diff << std::endl;

    _itr_count = i_itr;
    if (i_itr==max_itr) {
        std::cout << "the iteration reached the maximum count " << max_itr << std::endl;
    }

    return 0;
}

///// check if solution is converged
//template <class T_VARIABLE>
//bool AbstractIterativePartitionedMethod<T_VARIABLE>::isConverged(MyNamedVariableContainer& vars_prev, MyNamedVariableContainer& vars_current, double &v_diff)
//{
//    for (size_t i=0; i<vars_prev.size(); i++) {
//        double v_prev = .0;
//        //vars_prev.get(i)->eval(.0, v_prev);
//        double v_cur = .0;
//        //vars_current.get(i)->eval(.0, v_cur);
//        v_diff = std::abs(v_cur - v_prev);
//        if (v_diff>getEpsilon()) {
//            return false;
//        }
//    }
//    return true;
//}

}
