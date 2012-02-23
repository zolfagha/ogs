
#pragma once

#include <vector>
#include <map>
#include <algorithm>

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

    void create(size_t length, RowMajorSparsity *sparsity)
    {
        SparseTableCRS<IDX_TYPE> *crs = convertRowMajorSparsityToCRS<IDX_TYPE>(*sparsity);
        assert (length == crs->dimension);
        _A = new CRSMatrix<double, IDX_TYPE>(crs->dimension, crs->row_ptr, crs->col_idx, crs->data);
        _b.resize(length);
        _x.resize(length);
    }

    size_t getDimension() const
    {
        return _x.size();
    }

    virtual void reset()
    {
        (*_A) = .0;
        _b.assign(_b.size(), .0);
        _x.assign(_x.size(), .0);
        _vec_knownX_id.clear();
        _vec_knownX_x.clear();
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
                const size_t colId = vec_col_pos[j];
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

    void solve()
    {
        if (_vec_knownX_id.size()>0) {
            CRSMatrix<double, IDX_TYPE>* tmp_A = getA()->clone();
            double *org_eqsRHS = getRHS();
            double *org_eqsX = getX();
            std::vector<double> _tmp_b;
            std::vector<double> _tmp_x;
            std::map<size_t,size_t> _map_solved_orgEqs;

            //std::cout << "#before\n";
            //_tmp_A->printMat();
            setKnownXi_ReduceSizeOfEQS(tmp_A, org_eqsRHS, org_eqsX, _vec_knownX_id, _vec_knownX_x, _tmp_b, _tmp_x, _map_solved_orgEqs);
            //std::cout << "\n#after\n";
            //_tmp_A->printMat();

            solveEqs(tmp_A, &_tmp_b[0], &_tmp_x[0]);

            const size_t dim = tmp_A->getNRows();
            for (size_t i=0; i<dim; i++) {
                setX(_map_solved_orgEqs[i], _tmp_x[i]);
            }

            delete tmp_A;
        } else {
            solveEqs(getA(), getRHS(), getX());
        }
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

    void setKnownXi_ReduceSizeOfEQS(CRSMatrix<double, IDX_TYPE> *A, double *org_eqsRHS, double *org_eqsX, const std::vector<size_t> &vec_id, const std::vector<double> &vec_x, std::vector<double> &out_b, std::vector<double> &out_x, std::map<size_t,size_t> &map_solved_orgEqs)
    {
        assert(vec_id.size()==vec_x.size());

        const size_t n_org_rows = A->getNRows();

        std::vector<IDX_TYPE> removed_rows(vec_id.size());
        for (size_t i=0; i<vec_id.size(); i++) {
            const size_t id = vec_id[i];
            const double val = vec_x[i];
            removed_rows[i] = id;

            //b_i -= A(i,k)*val, i!=k
            for (size_t j=0; j<A->getNCols(); j++)
                org_eqsRHS[j] -= A->getValue(j, id)*val;
            //b_k = A_kk*val
            org_eqsRHS[id] = val; //=eqsA(id, id)*val;
            org_eqsX[id] = val; //=eqsA(id, id)*val;
        }

        //remove rows and columns
        std::sort(removed_rows.begin(), removed_rows.end());
        A->eraseEntries(removed_rows.size(), &removed_rows[0]);
        const size_t n_new_rows = n_org_rows-removed_rows.size();

        //remove X,RHS
        out_b.resize(n_new_rows);
        out_x.resize(n_new_rows);
        size_t new_id = 0;
        for (size_t i=0; i<n_org_rows; i++) {
            if (std::find(removed_rows.begin(), removed_rows.end(), i)!=removed_rows.end()) continue;
            out_b[new_id] = org_eqsRHS[i];
            out_x[new_id] = org_eqsX[i];
            map_solved_orgEqs[new_id] = i;
            new_id++;
        }
    }
};

/**
 * \brief 
 */
class SparseLinearEquations : public CRSLinearEquationsBase<unsigned>
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

    SparseLinearEquations() {};

    virtual ~SparseLinearEquations()
    {
    }

    void initialize() {};
    void finalize() {};

    void setOption(const Base::Options &option);

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

    DISALLOW_COPY_AND_ASSIGN(SparseLinearEquations);
};

}