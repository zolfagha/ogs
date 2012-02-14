
#pragma once

#include "NumLib/Discrete/DiscretizedEQS.h"

namespace NumLib
{

class CouplingEquations
{
private:

public:
    size_t addSystem(size_t var_id, std::vector<int> list_submatrix, std::vector<int> list_subRHS);

};

}
