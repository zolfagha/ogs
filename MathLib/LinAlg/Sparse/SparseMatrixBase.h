#ifndef SPARSEMATRIXBASE_H
#define SPARSEMATRIXBASE_H

#include "../MatrixBase.h"

namespace MathLib {

template<class T, class INT_TYPE> 
class TemplateSparseMatrixBase : public TemplateMatrixBase<INT_TYPE>
{
public:
	TemplateSparseMatrixBase(INT_TYPE n1, INT_TYPE n2) : TemplateMatrixBase<INT_TYPE> (n1,n2) {}
	TemplateSparseMatrixBase() : TemplateMatrixBase<INT_TYPE> () {}
	virtual void amux(T d, T const * const x, T *y) const = 0;         // y +=d*Ax
	virtual ~TemplateSparseMatrixBase() { }
};

template<class T> 
class SparseMatrixBase : public TemplateMatrixBase<unsigned>
{
public:
  SparseMatrixBase(unsigned n1, unsigned n2) : TemplateMatrixBase<INT_TYPE> (n1,n2) {}
  SparseMatrixBase() : TemplateMatrixBase<INT_TYPE> () {}
  virtual void amux(T d, T const * const x, T *y) const = 0;         // y +=d*Ax
  virtual ~SparseMatrixBase() { }
};


} // end namespace MathLib

#endif
