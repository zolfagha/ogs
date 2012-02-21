
#pragma once

#include "NumLib/Solution/TransientFEModel.h"
#include "NumLib/Solution/ISolution.h"

namespace ModelLib
{

#if 0
class GROUNDWATER_FLOW : public NumLib::TransientFemProblem
{
private:
    NumLib::SingleStepLinearFemAlgorithm _solution;
public:
    void initialiez();
    void setupModel(int*)
    {

    }

    bool solveTimeStep(NumLib::TimeStep t_n1)
    {
        return true;
    }

};
#endif


}
