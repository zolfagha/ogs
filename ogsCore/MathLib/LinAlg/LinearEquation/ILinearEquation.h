/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ILinearEquations.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

#include <iostream>
#include <vector>
#include "Eigen"
#include "BaseLib/Options.h"
#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

namespace MathLib
{

/**
 * \brief Interface for all algebraic linear equation systems
 */
class ILinearEquation
{
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;
    typedef Eigen::VectorXd LocalVector;

    virtual ~ILinearEquation() {};

    virtual void create(size_t length, RowMajorSparsity *sparsity=0) = 0;
    virtual bool isCreated() const = 0;
    virtual void setOption(const BaseLib::Options &option) = 0;
    virtual void reset() = 0;
    virtual size_t getDimension() const = 0;

    virtual double getA(size_t rowId, size_t colId) = 0;
    virtual void setA(size_t rowId, size_t colId, double v) = 0;
    virtual void addA(size_t rowId, size_t colId, double v) = 0;
    virtual void addAsub(const std::vector<size_t> &vec_row_pos, const std::vector<size_t> &vec_col_pos, LocalMatrix &sub_matrix, double fkt=1.0)
    {
        const size_t n_rows = vec_row_pos.size();
        const size_t n_cols = vec_col_pos.size();
        for (size_t i=0; i<n_rows; i++) {
            const size_t rowId = vec_row_pos[i];
            if (rowId==BaseLib::index_npos) continue;
            for (size_t j=0; j<n_cols; j++) {
                const size_t colId = vec_col_pos[j];
                if (colId==BaseLib::index_npos) continue;
                addA(rowId, colId, fkt*sub_matrix(i,j));
            }
        }
    }
    virtual void addAsub(std::vector<size_t> &vec_pos, LocalMatrix &sub_matrix, double fkt=1.0)
    {
        addAsub(vec_pos, vec_pos, sub_matrix, fkt);
    }

    virtual double getRHS(size_t rowId) = 0;
    virtual double* getRHS() = 0;
    virtual void setRHS(size_t rowId, double v) = 0;
    virtual void addRHS(size_t rowId, double v) = 0;
    virtual void addRHSsub(const std::vector<size_t> &vec_row_pos, double *sub_vector, double fkt=1.0)
    {
        for (size_t i=0; i<vec_row_pos.size(); i++) {
            const size_t rowId = vec_row_pos[i];
            if (rowId==BaseLib::index_npos) continue;
            addRHS(rowId, sub_vector[i]*fkt);
        }
    }

    virtual double* getX() = 0;

    virtual void setKnownX(size_t row_id, double x) = 0;
    virtual void setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x) = 0;
    virtual void solve() = 0;
    virtual void printout(std::ostream &os=std::cout) const = 0;
};

}
