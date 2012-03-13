
#pragma once

#include <vector>
#include <map>

#include "Base/CodingTools.h"

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
    CRSLinearEquationsBase() : _A(0) {};

    virtual ~CRSLinearEquationsBase()
    {
        Base::releaseObject(_A);
    }

    void create(size_t length, RowMajorSparsity *sparsity);
    bool isCreated() const { return _A!=0; };

    size_t getDimension() const
    {
        return _x.size();
    }

    virtual void reset();

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

    void addRHS(size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    double* getX()
    {
        return &_x[0];
    }

    void setX(size_t i, double v)
    {
        _x[i] = v;
    }

    void setKnownX(size_t id, double x)
    {
        _vec_knownX_id.push_back(id);
        _vec_knownX_x.push_back(x);
    }

    void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        _vec_knownX_id.insert(_vec_knownX_id.end(), vec_id.begin(), vec_id.end());
        _vec_knownX_x.insert(_vec_knownX_x.end(), vec_x.begin(), vec_x.end());
    }

    void solve();

    void printout(std::ostream &os=std::cout) const
    {
        os << "#A=" << std::endl;
        _A->printMat();
        os << "#x=" << std::endl;
        for (size_t i=0; i<_x.size(); i++) os << _x[i] << " "; os << std::endl;
        os << "#b=" << std::endl;
        for (size_t i=0; i<_b.size(); i++) os << _b[i] << " "; os << std::endl;
    }

protected:
    virtual void solveEqs(CRSMatrix<double, IDX_TYPE> *A, double *rhs, double *x) = 0;

private:
    CRSMatrix<double, IDX_TYPE> *_A;
    std::vector<double> _b;
    std::vector<double> _x;
    std::vector<size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;

    DISALLOW_COPY_AND_ASSIGN(CRSLinearEquationsBase);

    void setKnownXi_ReduceSizeOfEQS(CRSMatrix<double, IDX_TYPE> *A, double *org_eqsRHS, double *org_eqsX, const std::vector<size_t> &vec_id, const std::vector<double> &vec_x, std::vector<double> &out_b, std::vector<double> &out_x, std::map<size_t,size_t> &map_solved_orgEqs);
};


} //end
