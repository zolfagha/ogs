
#pragma once

#include <iostream>
#include <fstream>

namespace OGS5
{

class Matrix
{
public:
    Matrix(size_t rows, size_t cols = 1);
    Matrix();
    explicit Matrix(const Matrix& m);
    //
    void resize(size_t rows, size_t cols = 1);
    //
    virtual ~Matrix();
    void ReleaseMemory();                 //06.2010. WW

    // Operators
    virtual void operator= (double a);
    virtual void operator*= (double a);
    virtual void operator/= (double a);
    virtual void operator+= (double a);
    void operator= (const Matrix& m);
    void operator+= (const Matrix& m);
    void operator-= (const Matrix& m);

    void GetTranspose(Matrix& m);

    // vec_result = This*vec. vec_result must be initialized
    void multi(const double* vec, double* vec_result, double fac = 1.0);
    // m_result = this*m. m_result must be initialized
    void multi(const Matrix& m, Matrix& m_result, double fac = 1.0);
    // m_result = this*m1*m2. m_result must be initialized
    void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

    // Access to members
    virtual double& operator() (size_t i, size_t j = 0) const;
    void LimitSize(size_t nRows, size_t nCols = 1);

    size_t Rows() const {return nrows; }
    size_t Cols() const {return ncols; }
    size_t Size() const {return size; }

    // Print
    void Write(std::ostream& os = std::cout);
    void Write_BIN(std::fstream& os);
    void Read_BIN(std::fstream& is);
protected:
    size_t nrows, nrows0;
    size_t ncols, ncols0;
    size_t size;
    double* data;
    bool Sym;
};


} //end
