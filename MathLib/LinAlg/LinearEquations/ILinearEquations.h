
#pragma once

#include <vector>
#include "Base/Options.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

namespace MathLib
{

/**
 * \brief Interface for all algebraic linear equation systems
 */
class ILinearEquations
{
public:
    virtual void initialize() = 0;
    virtual void finalize() = 0;

    virtual void create(size_t length, RowMajorSparsity *sparsity=0) = 0;
    virtual void setOption(const Base::Options &option) = 0;
    virtual void reset() = 0;

    virtual double getA(size_t rowId, size_t colId) = 0;
    virtual void setA(size_t rowId, size_t colId, double v) = 0;
    virtual void addA(size_t rowId, size_t colId, double v) = 0;
    virtual void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0) = 0;
    virtual void addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0) = 0;

    virtual double getRHS(size_t rowId) = 0;
    virtual double* getRHS() = 0;
    virtual void setRHS(size_t rowId, double v) = 0;
    virtual void addRHS(std::vector<size_t> &vec_pos, double *sub_vector, double fkt=1.0) = 0;

    virtual double* getX() = 0;

    virtual void setKnownX(size_t row_id, double x) = 0;
    virtual void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x) = 0;
    virtual void solve() = 0;
};

}
