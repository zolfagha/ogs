
#pragma once

#include <vector>
#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"

namespace MathLib
{

/**
 * \brief 
 */
template<typename IDX_TYPE>
class CRSLinearEquationsBase : public ILinearEquations
{
public:
    void create(size_t length, RowMajorSparsity *sparsity)
    {
        SparseTableCRS<IDX_TYPE> *crs = convertRowMajorSparsityToCRS<IDX_TYPE>(*sparsity);
        assert (length == crs->dimension);
        _A = new CRSMatrix<double, IDX_TYPE>(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
        _b.resize(length);
        _x.resize(length);
    }

    void reset()
    {
        (*_A) = .0;
        _b.assign(_b.size(), .0);
        _x.assign(_x.size(), .0);
    }

    CRSMatrix<double, IDX_TYPE>* getA()
    {
        return _A;
    }

    double getA(size_t rowId, size_t colId)
    {
        return _A->getValue(rowId, colId);
    }

    void setA(size_t rowId, size_t colId, double v)
    {
        _A->setValue(rowId, colId, v);
    }

    void addA(size_t rowId, size_t colId, double v)
    {
        _A->addValue(rowId, colId, v);
    }

    void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    {
        for (size_t i=0; i<vec_row_pos.size(); i++) {
            const size_t rowId = vec_row_pos[i];
            for (size_t j=0; j<vec_col_pos.size(); j++) {
                const size_t colId = vec_row_pos[i];
                _A->addValue(rowId, colId, fkt*sub_matrix(i,j));
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

    }

    void setKnownX(std::vector<size_t> &vec_id, std::vector<double> &vec_x)
    {

    }


private:
    CRSMatrix<double, IDX_TYPE> *_A;
    std::vector<double> _b;
    std::vector<double> _x;
};

/**
 * \brief 
 */
class SparseLinearEquations : public CRSLinearEquationsBase<unsigned>
{
public:
    enum SolverType
    {
        CG,
        BiCGStab,
        GMRes
    };

    enum PreconditionerType
    {
        NONE,
        Diag,
        DiagScale
    };

    struct SpLinearOptions
    {
        SolverType solver_type;
        PreconditionerType precon_type;
        double error_tolerance;
        size_t max_iteration_step;

        SpLinearOptions()
        {
            solver_type = CG;
            precon_type = NONE;
            error_tolerance = 1.e-10;
            max_iteration_step = 500;
        };
    };

    void initialize() {};
    void finalize() {};

    //void create(size_t length, RowMajorSparsity *sparsity)
    //{
    //    SparseTableCRS<size_t> *crs = convertRowMajorSparsityToCRS<size_t>(*sparsity);
    //    assert (length == crs->dimension);
    //    _A = new CRSMatrix<double, size_t>(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
    //    _b.resize(length);
    //    _x.resize(length);

    //}

    void setOption(const Base::Options &option)
    {
        const Base::Options *op = option.getSubGroup("SpLinearOptions");
        if (op==0) return;

        if (op->hasOption("solver_type"))
            _option.solver_type = (SolverType)op->getOptionAsNum<int>("solver_type");
        if (op->hasOption("precon_type"))
            _option.precon_type = (PreconditionerType)op->getOptionAsNum<int>("precon_type");
        if (op->hasOption("error_tolerance"))
            _option.error_tolerance = op->getOptionAsNum<double>("error_tolerance");
        if (op->hasOption("max_iteration_step"))
            _option.max_iteration_step = op->getOptionAsNum<int>("max_iteration_step");
    }

    void setOption(const SpLinearOptions &option)
    {
        _option = option;
    }

    SpLinearOptions &getOption()
    {
        return _option;
    }

    //void reset()
    //{
    //    (*_A) = .0;
    //    _b.assign(_b.size(), .0);
    //    _x.assign(_x.size(), .0);
    //}

    //double getA(size_t rowId, size_t colId)
    //{
    //    return _A->getValue(rowId, colId);
    //}

    //void setA(size_t rowId, size_t colId, double v)
    //{
    //    _A->setValue(rowId, colId, v);
    //}

    //void addA(size_t rowId, size_t colId, double v)
    //{
    //    _A->addValue(rowId, colId, v);
    //}

    //void addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    //{
    //    for (size_t i=0; i<vec_row_pos.size(); i++) {
    //        const size_t rowId = vec_row_pos[i];
    //        for (size_t j=0; j<vec_col_pos.size(); j++) {
    //            const size_t colId = vec_row_pos[i];
    //            _A->addValue(rowId, colId, fkt*sub_matrix(i,j));
    //        }
    //    }
    //}

    //void addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt=1.0)
    //{
    //    addA(vec_pos, vec_pos, sub_matrix, fkt);
    //}

    //double getRHS(size_t rowId)
    //{
    //    return _b[rowId];
    //}

    //double* getRHS()
    //{
    //    return &_b[0];
    //}

    //void setRHS(size_t rowId, double v)
    //{
    //    _b[rowId] = v;
    //}

    //void addRHS(std::vector<size_t> &vec_row_pos, double *sub_vector, double fkt=1.0)
    //{
    //    for (size_t i=0; i<vec_row_pos.size(); i++) {
    //        const size_t rowId = vec_row_pos[i];
    //        _b[rowId] += sub_vector[i];
    //    }
    //}

    //double* getX()
    //{
    //    return &_x[0];
    //}

    //void setKnownX(size_t row_id, double x)
    //{

    //}

    //void setKnownX(std::vector<size_t> &vec_id, std::vector<double> &vec_x)
    //{

    //}

    void solve()
    {

    }

private:
    //CRSMatrix<double, size_t> *_A;
    //std::vector<double> _b;
    //std::vector<double> _x;
    SpLinearOptions _option;

};

/**
 * \brief 
 */
class DenseLinearEquations : public ILinearEquations
{
private:
    Matrix<double> *_A;
    std::vector<double> _b;
    std::vector<double> _x;

public:
    void initialize() {};
    void finalize() {};

    void create(size_t length, RowMajorSparsity *sparsity=0)
    {
        _A = new Matrix<double>(length, length);
        _b.resize(length);
        _x.resize(length);
    }

    void reset()
    {
        (*_A) = .0;
        _b.assign(_b.size(), .0);
        _x.assign(_x.size(), .0);
    }

    void setOption(const Base::Options &option)
    {

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

    }

    void setKnownX(std::vector<size_t> &vec_id, std::vector<double> &vec_x)
    {

    }

    void solve()
    {

    }

};

}