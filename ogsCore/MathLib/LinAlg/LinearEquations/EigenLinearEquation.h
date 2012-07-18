
#pragma once

#include <valarray>
#include <vector>
#include <Eigen>

#include "ILinearEquations.h"


namespace MathLib
{

/**
 * \brief Dense linear equation using Eigen
 */
class EigenDenseLinearEquation : public ILinearEquations
{
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
    typedef Eigen::VectorXd VectorType;

    EigenDenseLinearEquation() {};

    virtual ~EigenDenseLinearEquation() {};

    void create(size_t length, RowMajorSparsity* /*sp*/=0)
    {
        resize(length);
    }

    bool isCreated() const { return true; };

    void resize(size_t length)
    {
        _A.resize(length, length);
        _b.resize(length);
        _x.resize(length);
        reset();
    }

    void reset()
    {
        _A *= .0;
        _b *= .0;
        _x *= .0;
    }

    size_t getDimension() const { return _A.rows(); }

    MatrixType* getA() { return &_A; }

    double getA(size_t rowId, size_t colId)
    {
        return _A(rowId, colId);
    }

    void setA(size_t rowId, size_t colId, double v)
    {
        _A(rowId, colId) = v;
    }

    void addA(size_t rowId, size_t colId, double v)
    {
        _A(rowId, colId) += v;
    }

//    void addA(std::vector<size_t> &rowId, std::vector<size_t> &colId, MatrixType &m)
//    {
//        for (size_t i=0; i<rowId.size(); i++)
//            for (size_t j=0; j<colId.size(); j++)
//                _A(rowId[i], colId[j]) += m(i,j);
//    }
//
//    void addA(std::vector<size_t> &pos, MatrixType &m)
//    {
//        addA(pos, pos, m);
//    }

    double getRHS(size_t rowId)
    {
        return _b[rowId];
    }

    double* getRHS()
    {
        return &_b[0];
    }

    VectorType* getRHSAsVec() {return &_b;};

    void setRHS(size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    void addRHS(size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    double* getX()
    {
        return &_x[0];
    }

    VectorType* getXAsVec() {return &_x;};


    void setKnownX(size_t row_id, double x)
    {
        const size_t n_cols = _A.cols();
        //A(k, j) = 0.
        for (size_t j=0; j<n_cols; j++)
            _A(row_id, j) = .0;
        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<n_cols; j++)
            _b[j] -= _A(j, row_id)*x;
        //b_k = val
        _b[row_id] = x;
        //A(i, k) = 0., i!=k
        for (size_t j=0; j<n_cols; j++)
            _A(j, row_id) = .0;
        //A(k, k) = 1.0
        _A(row_id, row_id) = 1.0; //=x
    }

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        for (size_t i=0; i<vec_id.size(); ++i)
            setKnownX(vec_id[i], vec_x[i]);
    }

    void printout(std::ostream &os=std::cout) const
    {
        os << "not implemented yet." << std::endl;
    }

    void setOption(const BaseLib::Options &/*option*/) {};

    void solve()
    {
        _x = _A.colPivHouseholderQr().solve(_b);
    }

private:
    MatrixType _A;
    VectorType _b;
    VectorType _x;

};

}
