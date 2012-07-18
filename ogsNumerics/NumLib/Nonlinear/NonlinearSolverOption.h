
#pragma once

#include <string>

#include "BaseLib/CodingTools.h"

namespace NumLib
{

struct NonlinerSolverOption
{
    enum SolverType
    {
        Linear,
        Picard,
        Newton,
        INVALID
    };

    SolverType solver_type;
    double error_tolerance;
    long max_iteration;

    NonlinerSolverOption()
    {
        solver_type = Linear;
        error_tolerance = 1.e-6;
        max_iteration = 500;
    }

    SolverType getSolverType(const std::string &str)
    {
        RETURN_ENUM_IF_SAME_STRING(Linear, str);
        RETURN_ENUM_IF_SAME_STRING(Picard, str);
        RETURN_ENUM_IF_SAME_STRING(Newton, str);

        return INVALID;
    }
};

} //end

