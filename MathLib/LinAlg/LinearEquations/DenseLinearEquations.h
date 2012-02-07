
#pragma once

#include <vector>

#include "ILinearEquations.h"
#include "MathLib/LinAlg/Dense/Matrix.h"


namespace MathLib
{

/**
 * \brief 
 */
class DenseLinearEquationsBase : public ILinearEquations
{
public:
    typedef Matrix<double> MatrixType;
    typedef std::vector<double> VectorType;

    void create(size_t length, RowMajorSparsity *sparsity=0)
    {
        _A = new Matrix<double>(length, length);
        _b.resize(length);
        _x.resize(length);
        reset();
    }

    void reset()
    {
        (*_A) = .0;
        _b.assign(_b.size(), .0);
        _x.assign(_x.size(), .0);
    }

    size_t getDimension() const
    {
        return _A->getNRows();
    }

    MatrixType* getA() 
    {
        return _A;
    }

    double getA(size_t rowId, size_t colId)
    {
        return (*_A)(rowId, colId);
    }

    void setA(size_t rowId, size_t colId, double v)
    {
        (*_A)(rowId, colId) = v;
    }

    void addA(size_t rowId, size_t colId, double v)
    {
        (*_A)(rowId, colId) += v;
    }

    void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    {
        for (size_t i=0; i<vec_row_pos.size(); i++) {
            const size_t rowId = vec_row_pos[i];
            for (size_t j=0; j<vec_col_pos.size(); j++) {
                const size_t colId = vec_row_pos[i];
                (*_A)(rowId, colId) += fkt*sub_matrix(i,j);
            }
        }
    }

    void addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    {
        addA(vec_pos, vec_pos, sub_matrix, fkt);
    }

    double getRHS(size_t rowId)
    {
        return _b[rowId];
    }

    double* getRHS()
    {
        return &_b[0];
    }

    void setRHS(size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    void addRHS(std::vector<size_t> &vec_row_pos, double *sub_vector, double fkt=1.0)
    {
        for (size_t i=0; i<vec_row_pos.size(); i++) {
            const size_t rowId = vec_row_pos[i];
            _b[rowId] += sub_vector[i];
        }
    }

    double* getX()
    {
        return &_x[0];
    }

    void setKnownX(size_t row_id, double x)
    {
        const size_t n_cols = _A->getNCols();
        //A(k, j) = 0.
        for (size_t j=0; j<n_cols; j++)
            (*_A)(row_id, j) = .0;
        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<n_cols; j++)
            _b[j] -= (*_A)(j, row_id)*x;
        //b_k = val
        _b[row_id] = x;
        //A(i, k) = 0., i!=k
        for (size_t j=0; j<n_cols; j++)
            (*_A)(j, row_id) = .0;
        //A(k, k) = 1.0
        (*_A)(row_id, row_id) = 1.0; //=x
    }

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        for (size_t i=0; i<vec_id.size(); ++i)
            setKnownX(vec_id[i], vec_x[i]);
    }

private:
    Matrix<double> *_A;
    std::vector<double> _b;
    std::vector<double> _x;

};

/**
 * \brief 
 */
class DenseLinearEquations : public DenseLinearEquationsBase
{
public:
    void initialize() {};
    void finalize() {};

    void setOption(const Base::Options &option);

    void solve();
};


}