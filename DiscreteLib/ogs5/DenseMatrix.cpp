
#include "DenseMatrix.h"

namespace OGS5
{
// Constructors
Matrix::Matrix(size_t rows, size_t cols) :
	nrows (rows), nrows0 (rows), ncols (cols), ncols0 (cols),
	size (nrows * ncols), data (new double[size]), Sym (false)
{
	for(size_t i = 0; i < size; i++)
		data[i] = 0.0;
}

Matrix::Matrix() :
	nrows (0), nrows0 (0), ncols (0), ncols0 (0),
	size (nrows * ncols), data (NULL), Sym (false)
{}

Matrix::Matrix(const Matrix& m) :
	nrows (m.nrows), nrows0 (m.nrows), ncols (m.ncols), ncols0 (m.ncols),
	size (nrows * ncols), data (new double[size]), Sym (m.Sym)
{
	for(size_t i = 0; i < size; i++)
		data[i] = m.data[i];
}

void Matrix::resize(size_t rows, size_t cols)
{
	if (size > 0)
	{
		delete[] data;
		data = NULL;
	}

	Sym = false;
	nrows = rows;
	ncols = cols;
	nrows0 = rows;
	ncols0 = ncols;
	size = nrows * ncols;
	data = new double[size];
	for(size_t i = 0; i < size; i++)
		data[i] = 0.0;
}

Matrix::~Matrix()
{
	delete [] data;
	data = NULL;
}
// 06.2010. WW
void Matrix::ReleaseMemory()
{
	delete [] data;
	data = NULL;
}

void Matrix::operator= (double a)
{
	for(size_t i = 0; i < size; i++)
		data[i] = a;
}
void Matrix::operator *= (double a)
{
	for(size_t i = 0; i < size; i++)
		data[i] *= a;
}
void Matrix::operator /= (double a)
{
	for(size_t i = 0; i < size; i++)
		data[i] /= a;
}
void Matrix::operator += (double a)
{
	for(size_t i = 0; i < size; i++)
		data[i] += a;
}
//
void Matrix::operator = (const Matrix& m)
{
#ifdef gDEBUG
	if(nrows != m.Rows() || ncols != m.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
		abort();
	}
#endif
	for(size_t i = 0; i < nrows; i++)
		for(size_t j = 0; j < ncols; j++)
			data[i * ncols + j] = m(i,j);
}

//
void Matrix::operator += (const Matrix& m)
{
#ifdef gDEBUG
	if(nrows != m.Rows() || ncols != m.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
		abort();
	}
#endif
	for(size_t i = 0; i < nrows; i++)
		for(size_t j = 0; j < ncols; j++)
			data[i * ncols + j] += m(i,j);
}

//
void Matrix::operator -= (const Matrix& m)
{
#ifdef gDEBUG
	if(nrows != m.Rows() || ncols != m.Cols()) //Assertion, will be removed
	{
		std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
		abort();
	}
#endif
	for(size_t i = 0; i < nrows; i++)
		for(size_t j = 0; j < ncols; j++)
			data[i * ncols + j] -= m(i,j);
}
//
void Matrix::GetTranspose(Matrix& m)
{
#ifdef gDEBUG
	if(ncols != m.Rows() && nrows != m.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
		abort();
	}
#endif

	for(size_t i = 0; i < m.Rows(); i++)
		for(size_t j = 0; j < m.Cols(); j++)
			//          m(i,j) = data[j*ncols+i];
			m(i,j) = (*this)(j,i);
}
//
// m_results = this*m. m_results must be initialized
void Matrix::multi(const Matrix& m, Matrix& m_result, double fac)
{
#ifdef gDEBUG
	if(ncols != m.Rows() && nrows != m_result.Rows() && m.Cols() != m_result.Cols())
	{
		std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
		abort();
	}
#endif
	for(size_t i = 0; i < m_result.Rows(); i++)
		for(size_t j = 0; j < m_result.Cols(); j++)
		{
			if(Sym && (j > i))
				continue;
			// m_result(i,j) = 0.0;
			for(size_t k = 0; k < ncols; k++)
				//            m_result(i,j) += fac*data[i*ncols+k]*m(k,j);
				m_result(i,j) += fac * (*this)(i,k) * m(k,j);
		}
}

//
// m_results = this*m1*m2. m_results must be  initialized
void Matrix::multi(const Matrix& m1, const Matrix& m2, Matrix& m_result)
{
#ifdef gDEBUG
	if(ncols != m1.Rows() && m1.Cols() != m2.Rows()
	   && m2.Cols() != m_result.Cols() && nrows != m_result.Rows())
	{
		std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
		abort();
	}
#endif
	for(size_t i = 0; i < m_result.Rows(); i++)
		for(size_t j = 0; j < m_result.Cols(); j++)
		{
			if(Sym && (j > i))
				continue;
			//m_result(i,j) = 0.0;
			for(size_t k = 0; k < ncols; k++)
				for(size_t l = 0; l < m2.Rows(); l++)
					//                m_result(i,j) += data[i*ncols+k]*m1(k,l)*m2(l,j);
					m_result(i,j) += (*this)(i,k) * m1(k,l) * m2(l,j);
		}
}
// vec_result = This*vec. vec_result must be  initialized
void Matrix::multi(const double* vec, double* vec_result, double fac)
{
	for(int i = 0; (size_t)i < nrows; i++)
		for(int j = 0; (size_t)j < ncols; j++)
			vec_result[i] += fac * (*this)(i,j) * vec[j];
}

double& Matrix::operator() (size_t i, size_t j) const
{
#ifdef gDEBUG
	if(i >= nrows || j >= ncols)
	{
		std::cout << "\n Index exceeds the size of the matrix" << std::endl;
		abort();
	}
#endif
	return data[i * ncols + j];
}
void Matrix::LimitSize(size_t nRows, size_t nCols)
{
#ifdef gDEBUG
	if(nRows > nrows0 || nCols > ncols0)
	{
		std::cout << "\n Given size exceeds the original size of the matrix" << std::endl;
		abort();
	}
#endif
	nrows = nRows;
	ncols = nCols;
	size = nrows * ncols;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2004 WW Implementation
   02/2005 WW Change name
**************************************************************************/
void Matrix::Write(std::ostream& os)
{
	os.setf(std::ios::scientific, std::ios::floatfield);
	os.precision(12);

	for (size_t i = 0; i < nrows; i++)
	{
		os << "| ";
		for (size_t j = 0; j < ncols; j++)
			os << (*this)(i, j) << " ";
		os << "| " << std::endl;
	}
	os << std::endl;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   01/2006 WW Implementation
   03/2010 TF write whole matrix in one chunk
**************************************************************************/
void Matrix::Write_BIN(std::fstream& os)
{
	os.write((char*)data, size * sizeof(double));
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   01/2006 WW Implementation
**************************************************************************/
void Matrix::Read_BIN(std::fstream& is)
{
	for(size_t i = 0; i < size; i++)
		is.read((char*)(&data[i]), sizeof(data[i]));
}
} //end
