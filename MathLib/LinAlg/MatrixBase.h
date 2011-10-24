/*
 * MatrixBase.h
 *
 *  Created on: Sep 27, 2011
 *      Author: TF
 */

#ifndef MATRIXBASE_H_
#define MATRIXBASE_H_

namespace MathLib {

/**
 * class MatrixBase is the basis for all matrix classes (dense and sparse)
 */
template<typename INT_TYPE>
class TemplateMatrixBase {
public:
	/**
	 * Constructor for initialization of the number of rows and columns
	 * @param nrows number of rows
	 * @param ncols number of columns
	 * @return
	 */
	TemplateMatrixBase(INT_TYPE nrows=0, INT_TYPE ncols=0) :
		_n_rows(nrows), _n_cols(ncols)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
	TemplateMatrixBase (TemplateMatrixBase const& original) :
		_n_rows (original._n_rows), _n_cols (original._n_cols)
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~TemplateMatrixBase() {};
	/**
	 * get the number of rows
	 * @return the number of rows
	 */
	INT_TYPE getNRows () const { return _n_rows; }
	/**
	 * get the number of columns
	 * @return the number of columns
	 */
	INT_TYPE getNCols () const { return _n_cols; }

protected:
	/**
	 * the number of rows
	 */
	INT_TYPE _n_rows;
	/**
	 * the number of columns
	 */
	INT_TYPE _n_cols;
};

typedef TemplateMatrixBase<unsigned> MatrixBase;
}

#endif /* MATRIXBASE_H_ */
