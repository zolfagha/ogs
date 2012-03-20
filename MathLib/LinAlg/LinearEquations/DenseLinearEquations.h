
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

    DenseLinearEquationsBase() : _A(0) {};

    virtual ~DenseLinearEquationsBase()
    {
        if (_A) delete _A;
    }

    void create(size_t length, RowMajorSparsity* sp=0)
    {
        resize(length);
    }

    bool isCreated() const { return _A!=0; };


    void resize(size_t length);


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

    double getRHS(size_t rowId)
    {
        return _b[rowId];
    }

    double* getRHS()
    {
        return &_b[0];
    }

    VectorType* getRHSAsStdVec() {return &_b;};

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

    VectorType* getXAsStdVec() {return &_x;};


    void setKnownX(size_t row_id, double x);

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        for (size_t i=0; i<vec_id.size(); ++i)
            setKnownX(vec_id[i], vec_x[i]);
    }

    void printout(std::ostream &os=std::cout) const
    {
    	os << "not implemented yet." << std::endl;
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
