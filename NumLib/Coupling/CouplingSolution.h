
#pragma once

#include <cassert>
#include "SharedVariables.h"

namespace NumLib
{


/**
 * \brief Interface class of coupling problems
 */
class ICouplingProblem
{
public:
    virtual int solve() = 0;
    virtual void set(size_t varId, Variable*) = 0;
    virtual Variable* get(size_t varId) const = 0;
    ///// solve this problem
    //virtual int solve(SharedVariables* var=0) = 0;
    ///// update variables
    //virtual void update(SharedVariables* var) = 0;
    virtual bool check() = 0;

    virtual size_t getNumberOfInputVarameters() const = 0;
    virtual size_t getNumberOfOutputParameters() const = 0;
};

/**
 * \brief MonolithicSolution
 */
class IMonolithicProblem : public ICouplingProblem
{
public:
    void set(size_t in_var, Variable* var)
    {
        assert(in_var<_vec_in_var.size());
        _vec_in_var[in_var] = var;
    }

    void setInitial(size_t out, Variable* var)
    {
        assert(out<_vec_out_var.size());
        _vec_out_var[out] = var;
    }

    Variable* get(size_t out) const
    {
        assert(out<_vec_out_var.size());
        return _vec_out_var[out];
    }

    bool check() {return true;};

protected:
    std::vector<Variable*> _vec_in_var;
    std::vector<Variable*> _vec_out_var;
};

}
