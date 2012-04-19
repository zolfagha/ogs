
#pragma once

#include <string>
#include "Parameter.h"

namespace MathLib
{

/**
 * \brief Interface class having input and output parameters
 *
 * Concept
 * - users can set parameters as input or get parameters as output
 */
class IIOSystem
{
public:
    virtual ~IIOSystem() {};

    /// get the number of input parameters
    virtual size_t getNumberOfInputParameters() const = 0;

    /// set input parameter
    virtual void setInput(size_t id, const Parameter* val) = 0;

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameters() const = 0;

    /// get output parameter
    virtual const Parameter* getOutput(size_t id) const = 0;

protected:
    /// get input parameter
    virtual const Parameter* getInput(size_t id) const = 0;
};


} //end
