
#pragma once

#include "NamedVariableContainer.h"

namespace NumLib
{

/**
 * \brief Interface class of coupling problems
 */
class ICoupledProblem
{
public:
    virtual int solve() = 0;
    virtual void setParameter(size_t varId, Variable*) = 0;
    virtual Variable* getParameter(size_t varId) const = 0;
    virtual bool check() const= 0;
    virtual size_t getNumberOfInputParameters() const = 0;
    virtual size_t getNumberOfParameters() const = 0;
    virtual size_t getParameterIdForInput(size_t input_id) const = 0;
};

}
