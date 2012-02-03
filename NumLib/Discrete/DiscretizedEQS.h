
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"

namespace NumLib
{
//typedef MathLib::CRSMatrix<double,size_t> SparseMatrix;
//typedef MathLib::Matrix<double> DenseMatrix;
//typedef std::vector<double> Vectir;
//typedef MathLib::Matrix<double> LocalEQS;

/**
 * \brief Discrete equation
 */
typedef MathLib::SparseLinearEquations DiscretizedEQS;

#if 0
class DiscretizedEQS
{
private:
    MathLib::LinearEquations *ls;

public:

    SparseMatrix* getA() {return _A;};
    Vectir* getRHS() {return _RHS;}
    Vectir* getX() {return _X;};

    void create(size_t dimension)
    {
        throw std::exception("The method or operation is not implemented.");
    }

    void reset()
    {
        throw std::exception("The method or operation is not implemented.");
    }

    ///
    double getA(size_t row_id, size_t col_id);
    void setA(size_t row_id, size_t col_id, double v);
    void addA(size_t row_id, size_t col_id, double v);

    void add( std::vector<size_t> &dofmap, LocalEQS &localEQS ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }

    void setKnownX(size_t eqs_id, double x0);

    ///
    void solve();





};

#endif

}
