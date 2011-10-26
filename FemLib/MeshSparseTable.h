
#pragma once

#include <set>
#include "Mesh.h"
#include "LinAlg/Sparse/SparseTableCRS.h"

namespace FemLib
{

template<class INTTYPE>
MathLib::SparseTableCRS<INTTYPE>* generateSparseTableCRS(MeshLib::IMesh *mesh)
{
  const size_t n_nodes = mesh->getNumberOfNodes();

  //get number of nonzero 
  std::vector<INTTYPE> *ptr = new std::vector<INTTYPE>(n_nodes+1);
  std::vector<INTTYPE> *vec_col_idx = new std::vector<INTTYPE>();
  size_t counter_ptr = 0;
  size_t cnt_row = 0;
  for (size_t i=0; i<n_nodes; i++) {
    (*ptr)[cnt_row++] = counter_ptr;         // starting point of the row
    std::set<size_t> setConnection;
    MeshLib::Node *nod = mesh->getNode(i);
    setConnection.insert(i);
    // search connected nodes
    const std::set<size_t> &connected_nodes = nod->getConnectedNodes();
    for (std::set<size_t>::iterator it=connected_nodes.begin(); it!=connected_nodes.end(); it++) {
        setConnection.insert(*it);
    }
    for (std::set<size_t>::iterator it=setConnection.begin(); it!=setConnection.end(); it++) {
        vec_col_idx->push_back(*it);
        ++counter_ptr;
    }

  }
  (*ptr)[n_nodes] = counter_ptr;

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
  crs->row_ptr = &(*ptr)[0];
  crs->col_idx = &(*vec_col_idx)[0];
  crs->nonzero = vec_col_idx->size();
  crs->data = new double[vec_col_idx->size()];
  for (size_t i=0; i<vec_col_idx->size(); i++)
      crs->data[i] = .0;

  return crs;
}

}
