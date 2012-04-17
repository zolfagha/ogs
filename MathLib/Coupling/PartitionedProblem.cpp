
#include "PartitionedProblem.h"

namespace MathLib
{

int PartitionedProblem::solve()
{
    return _algorithm->solve(_list_subproblems, *getParameters(), _map);
}

} //end
