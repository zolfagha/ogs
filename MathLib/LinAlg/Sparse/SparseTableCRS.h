#pragma once

namespace MathLib {

template <class T>
struct SparseTableCRS { //temporally
    size_t dimension; 
    size_t nonzero;
    T *row_ptr;
    T *col_idx;
    double *data;
};

typedef struct SparseTableCRS<unsigned int> CRSUnsigned;
typedef struct SparseTableCRS<int> CRSSigned;

} // end namespace MathLib

