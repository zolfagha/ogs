#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <cassert>

//extern void CS_write(char*, unsigned, unsigned const*, unsigned const*, double const*);
//extern void CS_read(char*, unsigned&, unsigned*&, unsigned*&, double*&);

template<class T, class IDX_TYPE> void CS_write(std::ostream &os, IDX_TYPE n, IDX_TYPE const* iA, IDX_TYPE const* jA, T const* A)
{
	os.write((char*) &n, sizeof(IDX_TYPE));
	os.write((char*) iA, (n + 1) * sizeof(IDX_TYPE));
	os.write((char*) jA, iA[n] * sizeof(IDX_TYPE));
	os.write((char*) A, iA[n] * sizeof(T));
}

template<class T, class IDX_TYPE> void CS_read(std::istream &is, IDX_TYPE &n, IDX_TYPE* &iA, IDX_TYPE* &jA, T* &A)
{
	is.read((char*) &n, sizeof(IDX_TYPE));
	if (iA != NULL) {
		delete[] iA;
		delete[] jA;
		delete[] A;
	}
	iA = new IDX_TYPE[n + 1];
	assert(iA != NULL);
	is.read((char*) iA, (n + 1) * sizeof(IDX_TYPE));

	jA = new IDX_TYPE[iA[n]];
	assert(jA != NULL);
	is.read((char*) jA, iA[n] * sizeof(IDX_TYPE));

	A = new T[iA[n]];
	assert(A != NULL);
	is.read((char*) A, iA[n] * sizeof(T));

#ifndef NDEBUG
	// do simple checks
	if (iA[0] != 0) std::cerr << std::endl << "CRS matrix: array iA doesn't start with 0"
					<< std::endl;

	IDX_TYPE i = 0;
	while (i < iA[n] && jA[i] < n)
		++i;
	if (i < iA[n]) std::cerr << std::endl << "CRS matrix: the " << i
					<< "th entry of jA has the value " << jA[i] << ", which is out of bounds."
					<< std::endl;
#endif
}

#endif

