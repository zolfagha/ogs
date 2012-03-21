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

double normEuklid (double const * const vec, size_t n)
{
	return sqrt (scpr (vec, vec, n));
}

template<class T>
double norm1(T &v, size_t n)
{
	double s = .0;
	for (size_t i=0; i<n; i++)
		s += v[i];
	return sqrt(s);
};


} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
