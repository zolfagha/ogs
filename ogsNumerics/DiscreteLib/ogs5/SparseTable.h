
#pragma once

#include <iostream>

namespace OGS5
{
class CPARDomain;
/// Sparse matrix storage type //04.2011. WW
enum StorageType { CRS, JDS};

class SparseTable
{
public:
    //SparseTable(MeshLib::CFEMesh* a_mesh,
                //bool quadratic,
                //bool symm = false,
                //StorageType stype = JDS);
    SparseTable(CPARDomain &m_dom, bool quadratic, bool symm = false);
    ~SparseTable();
    void Write(std::ostream &os = std::cout);
private:
    bool symmetry;
    // Topology mapping from data array to matrix
    long* entry_column;
    long* num_column_entries;             // number of entries of each columns in sparse table
    long* row_index_mapping_n2o;          // Row index of sparse table to row index of matrix
    long* row_index_mapping_o2n;          // Inverse of last
    long* diag_entry;                     // Global index to the index of  entry_column
    long _size_entry_column;
    long max_columns;
    long _rows;
    StorageType storage_type; //04.2011. WW
    friend class CSparseMatrix;
};


}
