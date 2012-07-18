
#pragma once

#include <string>
#include "TXFunction.h"

namespace NumLib
{

class TXFunctionBuilder
{
public:
    ITXFunction* create(const std::string &name, double v)
    {
        if (name.compare("CONSTANT")==0) {
            return new TXFunctionConstant(v);
        }
        return NULL;
    }
};


}
