
#pragma once

#include "Mesh.h"
#include "LinAlg/Sparse/SparseTableCRS.h"

namespace FemLib
{

template<class INTTYPE>
MathLib::SparseTableCRS<INTTYPE>* generateSparseTableCRS(MeshLib::IMesh *mesh)
{
  const size_t n_nodes = mesh->getNumberOfNodes();

  //get number of nonzero 
  std::vector<INTTYPE> ptr(n_nodes+1);
  std::vector<INTTYPE> vec_col_idx;
  size_t counter_ptr = 0;
  size_t cnt_row = 0;
  for (size_t i=0; i<n_nodes; i++) {
    ptr[cnt_row++] = counter_ptr;         // starting point of the row
    // search connected nodes
    {
      size_t j_node = 0;
      vec_col_idx.push_back(j_node);
      ++counter_ptr;
    }
  }
  ptr[n_nodes] = counter_ptr;

#if 0
  //output CRS
  cout << "PTR:" << endl;
  for (size_t i=0; i<A.rows()+1; i++)
    cout << ptr[i] << ", "; 
  cout << endl;
  cout << "ColID:" << endl;
  for (size_t i=0; i<nonzero; i++)
    cout << col_idx[i] << ", "; 
  cout << endl;
  cout << "Data:" << endl;
  for (size_t i=0; i<nonzero; i++)
    cout << crs_data[i] << ", "; 
  cout << endl;
#endif

  MathLib::SparseTableCRS<INTTYPE> *crs(new MathLib::SparseTableCRS<INTTYPE>);
  crs->dimension = n_nodes;
  crs->row_ptr = &ptr[0];
  crs->col_idx = &vec_col_idx[0];
  crs->data = new double[vec_col_idx.size()];

  return crs;
}

}
