
#pragma once

#include "NumLib/Solution/TransientFEModel.h"
#include "NumLib/Solution/ISolution.h"

namespace ModelLib
{

class GROUNDWATER_FLOW : public NumLib::TransientFEModel
{
private:
    NumLib::TimeEulerSpFEMLinearSolution _solution;
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

}
