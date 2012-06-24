
#pragma once

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/LinearEquations/SparseLinearEquationBase.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"

namespace MathLib
{

/**
 * \brief 
 */
class SparseLinearEquation : public AbstractCRSLinearEquation<unsigned>
{
public:
    enum SolverType
    {
        SolverCG,
        SolverBiCGStab,
        SolverGMRes
    };

    enum PreconditionerType
    {
        NONE,
        PreconDiag,
        PreconDiagScale
    };

    struct SpLinearOptions
    {
        SolverType solver_type;
        PreconditionerType precon_type;
        double error_tolerance;
        size_t max_iteration_step;

        SpLinearOptions()
        {
            solver_type = SolverCG;
            precon_type = NONE;
            error_tolerance = 1.e-10;
            max_iteration_step = 500;
        };
    };

    SparseLinearEquation() {};

    virtual ~SparseLinearEquation()
    {
    }

    void initialize() {};
    void finalize() {};

    void setOption(const BaseLib::Options &option);

    void setOption(const SpLinearOptions &option)
    {
        _option = option;
    }

    SpLinearOptions &getOption()
    {
        return _option;
    }


protected:
    void solveEqs(CRSMatrix<double, unsigned> *A, double *rhs, double *x);

private:
    SpLinearOptions _option;

    DISALLOW_COPY_AND_ASSIGN(SparseLinearEquation);
};

}
