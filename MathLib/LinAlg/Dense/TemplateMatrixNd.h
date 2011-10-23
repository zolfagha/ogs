/*
 * \file Matrix.h
 *
 *  Created on: Mar 24, 2010
 *      Author: TF
 *  modified on:Jul 13, 2010
 *      HS & ZC
 */

#pragma once

#include <new>
#include <exception>
#include <stdexcept>
#include <iostream>

namespace MathLib {

template <class T, size_t NROWS, size_t NCOLS> 
class TemplateMatrixNd
{
public:
   TemplateMatrixNd () {};
   TemplateMatrixNd (const T& val);
   TemplateMatrixNd (const TemplateMatrixNd<T,NROWS,NCOLS> &src);

   ~TemplateMatrixNd () {};

   size_t getNRows () const { return NROWS; }
   size_t getNCols () const { return NCOLS; }

   inline T & operator() (size_t row, size_t col) throw (std::range_error);
   inline T & operator() (size_t row, size_t col) const throw (std::range_error);

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

   T& getData () { return data; }
   inline void fill(const T& val);

   virtual void operator= (const T& a);
   virtual void operator*= (const T& a);
   virtual void operator/= (const T& a);
   virtual void operator+= (const T& a);

private:
   // zero based addressing, but Fortran storage layout
   //inline size_t address(size_t i, size_t j) const { return j*rows+i; };
   // zero based addressing, C storage layout
   inline size_t address(size_t i, size_t j) const { return i*getNCols()+j; };

   T data[NROWS*NCOLS];
};

template<class T, size_t NROWS, size_t NCOLS> TemplateMatrixNd<T,NROWS,NCOLS>::TemplateMatrixNd (const T& val)
{
  (*this) = val;
};

template<class T, size_t NROWS, size_t NCOLS> T& TemplateMatrixNd<T,NROWS,NCOLS>::operator() (size_t row, size_t col)
throw (std::range_error)
{
  if ( (row >= getNRows()) | ( col >= getNCols()) )
    throw std::range_error ("Matrix: op() const range error");
  return data [address(row,col)];
}


template<class T, size_t NROWS, size_t NCOLS> T& TemplateMatrixNd<T,NROWS,NCOLS>::operator() (size_t row, size_t col) const
throw (std::range_error)
{
  if ( (row >= getNRows()) | ( col >= getNCols()) )
    throw std::range_error ("Matrix: op() const range error");
  return data [address(row,col)];
}


template<class T, size_t NROWS, size_t NCOLS> void TemplateMatrixNd<T,NROWS,NCOLS>::operator = (const T& val)
{
  for (size_t i=0; i<NROWS*NCOLS; i++)
    data[i] = val;
};

template<class T, size_t NROWS, size_t NCOLS> void TemplateMatrixNd<T,NROWS,NCOLS>::operator *= (const T& val)
{
  for (size_t i=0; i<NROWS*NCOLS; i++)
    data[i] *= val;
};

template<class T, size_t NROWS, size_t NCOLS> void TemplateMatrixNd<T,NROWS,NCOLS>::operator /= (const T& val)
{
  for (size_t i=0; i<NROWS*NCOLS; i++)
    data[i] /= val;
};

template<class T, size_t NROWS, size_t NCOLS> void TemplateMatrixNd<T,NROWS,NCOLS>::operator += (const T& val)
{
  for (size_t i=0; i<NROWS*NCOLS; i++)
    data[i] += val;
};

typedef TemplateMatrixNd<double,3,3> Matrix3d;

}
