
#pragma once

#include <vector>
#include <iostream>

#include "ICoupledProblem.h"


namespace NumLib
{
typedef std::pair<ICoupledProblem*,size_t> PairSysVarId;
typedef std::pair<size_t, size_t> PairInputVar;

/**
 * \brief Interface class for partitioned algorithm
 */
class IPartitionedAlgorithm
{
public:
    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<ICoupledProblem*> &subproblems, NamedVariableContainer &vars, VariableMappingTable &mapping) = 0;
};

/**
 * \brief Abstract class for iterative partitioned methods
 */
class AbstractIterativePartitionedMethod : public IPartitionedAlgorithm
{
public:
    AbstractIterativePartitionedMethod() : _max_itr(100), _epsilon(1e-3)
    {
    }

    AbstractIterativePartitionedMethod(double epsilon, size_t max_count) : _max_itr(max_count), _epsilon(epsilon)
    {
    }

    /// get max. iteration
    size_t getMaximumIterationCounts() const {return _max_itr;};
    /// get epsilon
    double getEpsilon() const {return _epsilon;};
    /// get iteration count
    size_t getIterationCounts() const {return _itr_count;};

    /// solve
    int solve(std::vector<ICoupledProblem*> &subproblems, NamedVariableContainer &vars, VariableMappingTable &mapping)
    {
        const size_t n_subproblems = subproblems.size();

        // initialize variables
        const size_t n_para = vars.size();
        for (size_t i=0; i<n_para; i++) {
            PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            ICoupledProblem *solution = sysvarid.first;
            size_t internalId = sysvarid.second;
            if (solution!=0)
                vars.set(i, *const_cast<Variable*>(solution->getParameter(internalId)));
        }

        // obj to store previous state
        NamedVariableContainer vars_prev;

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
                ICoupledProblem *solution = subproblems[i];
                std::vector<PairInputVar> &list_input_var = mapping._list_subproblem_input_source[i];
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
                is_converged = isConverged(vars_prev, vars, v_diff);
                //std::cout << i_itr << ": diff=" << v_diff << std::endl;
            } else {
                is_converged = true;
            }
            ++i_itr;
        } while (!is_converged && i_itr<max_itr);
        std::cout << "iteration count=" << i_itr << ", error=" << v_diff << std::endl;

        _itr_count = i_itr;
        if (i_itr==max_itr) {
            std::cout << "the iteration reached the maximum count " << max_itr << std::endl;
        }

        return 0;
    }

protected:
    /// check if solution is converged
    bool isConverged(NamedVariableContainer& vars_prev, NamedVariableContainer& vars_current, double &v_diff)
    {
        for (size_t i=0; i<vars_prev.size(); i++) {
            double v_prev = vars_prev.get(i)->eval(.0);
            double v_cur = vars_current.get(i)->eval(.0);
            v_diff = fabs(v_cur - v_prev);
            if (v_diff>getEpsilon()) {
                return false;
            }
        }
        return true;
    }


    virtual void doPostAfterSolve( ICoupledProblem& solution, NamedVariableContainer& vars, VariableMappingTable &mapping )  {}

    virtual void doPostAfterSolveAll( NamedVariableContainer &vars, VariableMappingTable &mapping ) {}



private:
    size_t _max_itr;
    size_t _itr_count;
    double _epsilon;
};


/**
 * \brief Block Jacobi iterative partitioned method
 */
class BlockJacobiMethod : public AbstractIterativePartitionedMethod
{
public:
    BlockJacobiMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(epsilon, max_count)
    {
    }

    virtual void doPostAfterSolveAll( NamedVariableContainer &vars, VariableMappingTable &mapping ) 
    {
        // set current state to shared variables
        const size_t n_vars = vars.size();
        for (size_t i=0; i<n_vars; i++) {
            VariableMappingTable::PairSysVarId &v = mapping._map_paraId2subproblem[i];
            const ICoupledProblem *solution = v.first;
            if (solution!=0)
                vars.set(i, *const_cast<Variable*>(solution->getParameter(v.second)));
        }
    }

};


/**
 * \brief Block Gauss-Seidel iterative partitioned method
 */
class BlockGaussSeidelMethod : public AbstractIterativePartitionedMethod
{
public:
    BlockGaussSeidelMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(epsilon, max_count)
    {
    }

    virtual void doPostAfterSolve( ICoupledProblem & solution, NamedVariableContainer& vars, VariableMappingTable &mapping ) 
    {
        // update shared variables
        const size_t n_vars = vars.size();
        for (size_t i=0; i<n_vars; i++) {
            VariableMappingTable::PairSysVarId &v = mapping._map_paraId2subproblem[i];
            const ICoupledProblem *work_solution = v.first;
            if (work_solution==&solution) {
                vars.set(i, *const_cast<Variable*>(solution.getParameter(v.second)));
            }
        }
    }
};

/**
 * \brief 
 */
class ITransientPartitionedAlgorithm
{
public:
    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<ICoupledProblem*> &subproblems, NamedVariableContainer &vars_t_n, NamedVariableContainer &vars_t_n1, VariableMappingTable &mapping) = 0;
};

/**
 * \brief 
 */
template <class T_PARTITIONED>
class AbstractPartitionedStaggeredMethod : public ITransientPartitionedAlgorithm
{
public:
    AbstractPartitionedStaggeredMethod() : _part_method() {};
    AbstractPartitionedStaggeredMethod(double epsilon, size_t max_count) : _part_method(epsilon, max_count) {};

    /// solve
    int solve(std::vector<ICoupledProblem*> &subproblems, NamedVariableContainer &vars_t_n, NamedVariableContainer &vars_t_n1, VariableMappingTable &mapping)
    {
        // initialize variables at t_n1 from t_n
        const size_t n_para = vars_t_n.size();
        for (size_t i=0; i<n_para; i++) {
            PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            ICoupledProblem *solution = sysvarid.first;
            size_t internalId = sysvarid.second;
            if (solution!=0)
                vars_t_n1.set(i, *vars_t_n.get(i));
        }

        _part_method.solve(subproblems, vars_t_n1, mapping);

        return 0;
    }
private:
    T_PARTITIONED _part_method;
};

typedef AbstractPartitionedStaggeredMethod<BlockJacobiMethod> ParallelStaggeredMethod;
typedef AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod> SerialStaggeredMethod;

}
