
#pragma once

#include "ParameterTable.h"

namespace MathLib
{

/**
 * \brief Interface class of coupling problems
 */
class ICoupledSystem
{
public:

	virtual ~ICoupledSystem() {};

    /// solve
    virtual int solve() = 0;

    /// check consistency
    virtual bool check() const= 0;

    // parameter for both input and output
    virtual size_t getNumberOfParameters() const = 0;
    virtual const Parameter* getOutput(size_t parameter_id) const = 0;
    //virtual void setParameter(size_t parameter_id, Parameter*) = 0;

    // input parameters
    virtual size_t getNumberOfInputParameters() const = 0;
    virtual size_t getParameterIdForInput(size_t input_id) const = 0;
    virtual void setInput(size_t parameter_id, const Parameter* val) = 0;

    // output parameters
    //virtual size_t getNumberOfOutputParameters() const = 0;
    //virtual size_t getParameterIdForOutput(size_t input_id) const = 0;

    //virtual void setOutput(size_t parameter_id, Parameter* val) = 0;
protected:
    virtual const Parameter* getInput(size_t parameter_id) const = 0;
};

}
