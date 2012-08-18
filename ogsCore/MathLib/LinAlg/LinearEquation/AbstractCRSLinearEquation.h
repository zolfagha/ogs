/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractCRSLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "ILinearEquation.h"

namespace MathLib
{
/**
 * \brief 
 */
template<typename IDX_TYPE>
class AbstractCRSLinearEquation : public ILinearEquation
{
public:
    typedef CRSMatrix<double, IDX_TYPE> MatrixType;
    typedef std::vector<double> VectorType;
    
    AbstractCRSLinearEquation() : _A(0) {};

    virtual ~AbstractCRSLinearEquation()
    {
        BaseLib::releaseObject(_A);
    }

    void create(size_t length, RowMajorSparsity *sparsity);
    bool isCreated() const { return _A!=0; };

    size_t getDimension() const
    {
        return _x.size();
    }

    virtual void reset();

    MatrixType* getA()
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
    virtual void solveEqs(MatrixType *A, double *rhs, double *x) = 0;

private:
    MatrixType *_A;
    VectorType _b;
    VectorType _x;
    std::vector<size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;

    DISALLOW_COPY_AND_ASSIGN(AbstractCRSLinearEquation);

    void setKnownXi_ReduceSizeOfEQS(MatrixType *A, double *org_eqsRHS, double *org_eqsX, const std::vector<size_t> &vec_id, const std::vector<double> &vec_x, std::vector<double> &out_b, std::vector<double> &out_x, std::map<size_t,size_t> &map_solved_orgEqs);
};


} //end
