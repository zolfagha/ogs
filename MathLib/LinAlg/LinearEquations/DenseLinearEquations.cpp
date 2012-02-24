
#include "DenseLinearEquations.h"

#include <algorithm>

#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"


namespace MathLib
{

void DenseLinearEquationsBase::resize(size_t length)
{
    if (_A==0) {
        _A = new Matrix<double>(length, length);
    } else {
        _A->resize(length, length);
    }
    _b.resize(length);
    _x.resize(length);
    reset();
}

void DenseLinearEquationsBase::setKnownX(size_t row_id, double x)
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

void DenseLinearEquations::setOption(const Base::Options &option)
{

}

void DenseLinearEquations::solve()
{
    MathLib::GaussAlgorithm solver(*this->getA());
    DenseLinearEquationsBase::VectorType *b = this->getRHSAsStdVec();
    DenseLinearEquationsBase::VectorType *x = this->getXAsStdVec();
    std::copy(b->begin(), b->end(), x->begin());
    solver.execute(&(*x)[0]);
}


} //end namespace
