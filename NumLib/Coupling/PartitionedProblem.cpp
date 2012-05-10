
#include "PartitionedProblem.h"

namespace NumLib
{

int PartitionedProblem::solve()
{
    return _algorithm->solve(_list_subproblems, *getParameters(), _map);
}

} //end
