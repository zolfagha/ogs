/*
 * VectorNorms.h
 *
 *  Created on: Jun 6, 2011
 *      Author: TF
 */

#ifndef VECTORNORMS_H_
#define VECTORNORMS_H_

#include <cmath>

#include "MathLib/MathTools.h"

namespace MathLib {

inline double normEuklid (double const * const vec, size_t n)
{
	return sqrt (scpr (vec, vec, n));
}

template<class T>
inline double norm_p(T &v, size_t n, int p)
{
	double s = .0;
	for (size_t i=0; i<n; i++)
		s += pow(fabs(v[i]),p);
	s = pow(s, 1./p);
	return s;
};

template<class T>
inline double norm1(T &v, size_t n)
{
	double s = .0;
	for (size_t i=0; i<n; i++)
		s += fabs(v[i]);
	return s;
};


} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
