/*
 * CRSMatrix.h
 *
 *  Created on: Sep 20, 2011
 *      Author: TF
 */

#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <string>
#include <fstream>
#include <iostream>
#include "SparseMatrixBase.h"
#include "sparse.h"
#include "amuxCRS.h"
#include "../Preconditioner/generateDiagPrecond.h"
#include "binarySearch.h"

namespace MathLib {

template<class T, class INT_TYPE> 
class TemplateCRSMatrix : public TemplateSparseMatrixBase<T, INT_TYPE>
{
public:
	TemplateCRSMatrix(std::string const &fname) :
		TemplateSparseMatrixBase<T, INT_TYPE>(),
		_number_of_nonzero(0), _row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{
		std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
		if (in) {
			CS_read(in, TemplateSparseMatrixBase<T, INT_TYPE>::_n_rows, _row_ptr, _col_idx, _data);
			TemplateSparseMatrixBase<T, INT_TYPE>::_n_cols = TemplateSparseMatrixBase<T, INT_TYPE>::_n_rows;
			in.close();
		} else {
			std::cout << "cannot open " << fname << std::endl;
		}
	}

	TemplateCRSMatrix(INT_TYPE n, INT_TYPE *iA, INT_TYPE *jA, T* A) :
		TemplateSparseMatrixBase<T, INT_TYPE>(n,n), _number_of_nonzero(0), _row_ptr(iA), _col_idx(jA), _data(A)
	{}

    TemplateCRSMatrix(INT_TYPE n, INT_TYPE number_of_nonzero, INT_TYPE *iA, INT_TYPE *jA, T* A) :
    TemplateSparseMatrixBase<T, INT_TYPE>(n,n), _number_of_nonzero(number_of_nonzero), _row_ptr(iA), _col_idx(jA), _data(A)
    {}
        
	TemplateCRSMatrix(INT_TYPE n1) :
		TemplateSparseMatrixBase<T, INT_TYPE>(n1, n1), _row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{}

	virtual ~TemplateCRSMatrix()
	{
		delete [] _row_ptr;
		delete [] _col_idx;
		delete [] _data;
	}

	virtual void amux(T d, T const * const x, T *y) const
	{
      std::cout << "Not implemented." << std::endl;
		//amuxCRS(d, TemplateMatrixBase::_n_rows, _row_ptr, _col_idx, _data, x, y);
	}

    virtual void precondApply(T* x) const
    {}

    virtual INT_TYPE getNumberOfNonzeroEntries() const 
    {
       return _number_of_nonzero;
    }

    virtual INT_TYPE* getRowPtr() {
        return _row_ptr;
    } 

    virtual INT_TYPE* getColIdx() {
        return _col_idx;
    } 

    virtual T* getRawData() {
        return _data;
    } 

    inline T & operator() (size_t row, size_t col)
    {
#if 0
      if ( (row >= getNRows()) | ( col >= getNCols()) )
        throw std::range_error ("Matrix: op() const range error");
#endif
      const size_t _row_i_begin = _row_ptr[row];
      setZero();
      if (_row_i_begin>=_number_of_nonzero)
        return zero_e;

      size_t col_id_begin = _col_idx[_row_i_begin];
      size_t col_id_end = 0;
      if (row==this->getNRows()-1)
        col_id_end = _col_idx[_number_of_nonzero-1];
      else
        col_id_end = _col_idx[_row_ptr[row+1]];

      if (col < col_id_begin || col_id_end < col)
        return zero_e;

      size_t k = searchElement<INT_TYPE, T*>(col, col_id_begin, col_id_end, _data); 
      if(k==-1) {
        return zero_e; 
      }

      return _data[k]; // 
    }


protected:
	INT_TYPE *_row_ptr;
	INT_TYPE *_col_idx;
    INT_TYPE _number_of_nonzero;
	T* _data;
    mutable T zero_e; 
    inline void setZero();
};

template<class T, class INT_TYPE> void TemplateCRSMatrix<T, INT_TYPE>::setZero()
{
  zero_e = .0;
};

//template<class INT_TYPE> void TemplateCRSMatrix<double, INT_TYPE>::setZero()
//{
//  zero_e = .0;
//};

template<class T>
class CRSMatrix : public TemplateCRSMatrix<T, unsigned>
{
public:
	CRSMatrix(std::string const &fname) :
		TemplateCRSMatrix<T, unsigned>(fname)
	{
	}

	CRSMatrix(unsigned n, unsigned *iA, unsigned *jA, T* A) :
		TemplateCRSMatrix<T, unsigned> (n, iA, jA, A)
	{}

    CRSMatrix(unsigned n, unsigned number_of_nonzero, unsigned *iA, unsigned *jA, T* A) :
        TemplateCRSMatrix<T, unsigned> (n, number_of_nonzero, iA, jA, A)
    {}
        
	CRSMatrix(unsigned n1) :
        TemplateCRSMatrix<T, unsigned> (n1)
	{}

	virtual ~CRSMatrix()
	{
	}

    virtual void amux(T d, T const * const x, T *y) const
    {
      amuxCRS(d, this->getNRows(), TemplateCRSMatrix<T, unsigned>::_row_ptr, TemplateCRSMatrix<T, unsigned>::_col_idx, TemplateCRSMatrix<T, unsigned>::_data, x, y);
//      amuxCRS(d, TemplateCRSMatrix<T, unsigned>::_n_rows, _row_ptr, _col_idx, _data, x, y);
    }
};

} // end namespace MathLib

#endif

