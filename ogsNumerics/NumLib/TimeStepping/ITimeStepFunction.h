
#pragma once

#include "TimeUnit.h"

namespace NumLib
{

class ITimeStepFunction
{
public:
    virtual double getBeginning() const = 0;
    virtual double getEnd() const = 0;
    virtual double getPrevious() const = 0;
    virtual double getNext(double t_current) = 0;
    virtual void accept() = 0;
    virtual ITimeStepFunction* clone() = 0;
    virtual ~ITimeStepFunction() {};
};

}
