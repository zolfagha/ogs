
#include "SparseTable.h"

#include "par_ddc_group.h"

namespace OGS5
{
#if 0
/*\!
 ********************************************************************
   Create sparse matrix table
   01/2006 WW
   08/2007 WW
   10/2007 WW
   02/2008 PCH  Compressed Row Storage with LIS option
   03/2010 WW: CRS storage for matrix algebraic
 ********************************************************************
 */
SparseTable::SparseTable(CFEMesh* a_mesh, bool quadratic, bool symm, StorageType stype)
    : symmetry(symm), storage_type(stype)
{
    long i = 0, j = 0, ii = 0, jj = 0;
    long lbuff0 = 0, lbuff1 = 0;
    long** larraybuffer;
    larraybuffer = NULL;
    //
    // In sparse table, = number of nodes
    rows = a_mesh->GetNodesNumber(quadratic);
    size_entry_column = 0;
    diag_entry = new long[rows];

    if(storage_type == JDS)
    {
        row_index_mapping_n2o = new long[rows];
        row_index_mapping_o2n = new long[rows];
    }
    else if (storage_type == CRS)
    {
        row_index_mapping_n2o = NULL;
        row_index_mapping_o2n = NULL;
    }

    if(symmetry)
    {
        larraybuffer = new long*[rows];
        for(i = 0; i < rows; i++)
        {
            if(storage_type == JDS)
                row_index_mapping_n2o[i] = i;
            // 'diag_entry' used as a temporary array
            // to store the number of nodes connected to this node
            lbuff1 = (long)a_mesh->nod_vector[i]->getConnectedNodes().size();
            larraybuffer[i] = new long[lbuff1 + 1];
            //
            larraybuffer[i][0] = lbuff1;
            for(j = 0; j < lbuff1; j++)
                larraybuffer[i][j +
                                1] = a_mesh->nod_vector[i]->getConnectedNodes()[j];
            a_mesh->nod_vector[i]->getConnectedNodes().clear();
            for(j = 0; j < lbuff1; j++)
            {
                jj = larraybuffer[i][j + 1];
                if(i <= jj)
                    a_mesh->nod_vector[i]->getConnectedNodes().push_back(jj);
            }
        }
    }

    /// CRS storage
    if(storage_type == CRS)
    {
        /// num_column_entries saves vector ptr of CRS
        num_column_entries = new long[rows + 1];

        std::vector<long> A_index;
        long col_index;

        for(i = 0; i < rows; i++)
        {
            num_column_entries[i] = (long)A_index.size();

            for(j = 0; j < (long)a_mesh->nod_vector[i]->getConnectedNodes().size(); j++)
            {
                col_index = a_mesh->nod_vector[i]->getConnectedNodes()[j];

                /// If linear element is used
                if((!quadratic) && (col_index >= rows))
                    continue;

                if(i == col_index)
                    diag_entry[i] = (long)A_index.size();
                A_index.push_back(col_index);
            }
        }

        size_entry_column = (long)A_index.size();
        num_column_entries[rows] = size_entry_column;

        entry_column = new long[size_entry_column];
        for(i = 0; i < size_entry_column; i++)
            entry_column[i] = A_index[i];
    }
    else if(storage_type == JDS)
    {
        //
        //--- Sort, from that has maximum connect nodes to that has minimum connect nodes
        //
        for(i = 0; i < rows; i++)
        {
            row_index_mapping_n2o[i] = i;
            // 'diag_entry' used as a temporary array
            // to store the number of nodes connected to this node
            diag_entry[i] = (long)a_mesh->nod_vector[i]->getConnectedNodes().size();
            if(!quadratic)
            {
                lbuff0 = 0;
                for(j = 0; j < diag_entry[i]; j++)
                    if(a_mesh->nod_vector[i]->getConnectedNodes()[j] <
                       static_cast<size_t>(rows))
                        lbuff0++;
                diag_entry[i] = lbuff0;
            }
            size_entry_column += diag_entry[i];
        }

        //
        for(i = 0; i < rows; i++)
        {
            // 'diag_entry' used as a temporary array
            // to store the number of nodes connected to this node
            lbuff0 = diag_entry[i]; // Nodes to this row
            lbuff1 = row_index_mapping_n2o[i];
            j = i;
            while((j > 0) && (diag_entry[j - 1] < lbuff0))
            {
                diag_entry[j] = diag_entry[j - 1];
                row_index_mapping_n2o[j] = row_index_mapping_n2o[j - 1];
                j = j - 1;
            }
            diag_entry[j] = lbuff0;
            row_index_mapping_n2o[j] = lbuff1;
        }
        // Old index to new one
        for(i = 0; i < rows; i++)
            row_index_mapping_o2n[row_index_mapping_n2o[i]] = i;
        // Maximum number of columns in the sparse table
        max_columns = diag_entry[0];
        //--- End of sorting
        //
        //--- Create sparse table
        //
        num_column_entries = new long[max_columns];
        entry_column = new long[size_entry_column];
        // 1. Count entries in each column in sparse table
        for (i = 0; i < max_columns; i++)
            num_column_entries[i] = 0;
        for (i = 0; i < rows; i++)
            // 'diag_entry' still is used as a temporary array
            // it stores that numbers of nodes connect to this nodes
            for (j = 0; j < diag_entry[i]; j++)
                num_column_entries[j]++;

        // 2. Fill the sparse table, i.e. store all its entries to
        //    entry_column
        lbuff0 = 0;

        for (i = 0; i < max_columns; i++)
            for (j = 0; j < num_column_entries[i]; j++)
            {
                // ii is the real row index of this entry in matrix
                ii = row_index_mapping_n2o[j];
                // jj is the real column index of this entry in matrix
                jj = a_mesh->nod_vector[ii]->getConnectedNodes()[i];
                entry_column[lbuff0] = jj;

                // Till to this stage, 'diag_entry' is really used to store indices of the diagonal entries.
                // Hereby, 'index' refers to the index in entry_column array.
                if(ii == jj)
                    diag_entry[ii] = lbuff0;
                //
                lbuff0++;
            }
    }

    // For the case of symmetry matrix
    if(symmetry)
    {
        for(i = 0; i < rows; i++)
        {
            lbuff0 = larraybuffer[i][0];
            a_mesh->nod_vector[i]->getConnectedNodes().resize(lbuff0);
            //
            for(j = 0; j < lbuff0; j++)
                a_mesh->nod_vector[i]->getConnectedNodes()[j] =
                        larraybuffer[i][j + 1];
        }
        for(i = 0; i < rows; i++)
        {
            delete [] larraybuffer[i];
            larraybuffer[i] = 0;
        }
        delete [] larraybuffer;
        larraybuffer = 0;
    }
}
#endif

/*\!
 ********************************************************************
   Create sparse matrix table for each domain
   12/2007 WW
 ********************************************************************
 */
SparseTable::SparseTable(CPARDomain &m_dom, bool quadratic, bool symm) :
    symmetry(symm)
{
    long i = 0, j = 0, ii = 0, jj = 0;
    long lbuff0 = 0, lbuff1 = 0;
    storage_type = JDS;
    //
    _rows = m_dom.getNumberOfDomainNodes(quadratic);
    _size_entry_column = 0;
    //
    row_index_mapping_n2o = new long[_rows];
    row_index_mapping_o2n = new long[_rows];
    diag_entry = new long[_rows];

    if(symmetry)
    {
        std::vector<long> conc;

        for(i = 0; i < _rows; i++)
        {
            row_index_mapping_n2o[i] = i;
            // 'diag_entry' used as a temporary array
            // to store the number of nodes connected to this node
            lbuff1 = m_dom.getNumberOfNodesConnectedToNode(i);
            //
            for(j = 0; j < lbuff1; j++)
            {
                jj = m_dom.get_node_conneted_nodes(i,j);
                if(i <= jj)
                    conc.push_back(jj);
                m_dom.set_node_conneted_nodes(i,j,0);
            }
            // Number of nodes connected to this node.
            m_dom.setNumberOfNodesConnectedToNode(i, (long)conc.size());
            // New
            for(j = 0; j < m_dom.getNumberOfNodesConnectedToNode(i); j++)
                m_dom.set_node_conneted_nodes(i, j, conc[j]);
        }
    }
    //
    //--- Sort, from that has maximum connect nodes to that has minimum connect nodes
    //
    for(i = 0; i < _rows; i++)
    {
        row_index_mapping_n2o[i] = i;
        // 'diag_entry' used as a temporary array
        // to store the number of nodes connected to this node
        diag_entry[i] =  m_dom.getNumberOfNodesConnectedToNode(i);
        if(!quadratic)
        {
            lbuff0 = 0;
            for(j = 0; j < diag_entry[i]; j++)
                if(m_dom.get_node_conneted_nodes(i, j) < _rows)
                    lbuff0++;
            diag_entry[i] = lbuff0;
        }
        _size_entry_column += diag_entry[i];
    }
    //
    for(i = 0; i < _rows; i++)
    {
        // 'diag_entry' used as a temporary array
        // to store the number of nodes connected to this node
        lbuff0 = diag_entry[i];   // Nodes to this row
        lbuff1 = row_index_mapping_n2o[i];
        j = i;
        while((j > 0) && (diag_entry[j - 1] < lbuff0))
        {
            diag_entry[j] = diag_entry[j - 1];
            row_index_mapping_n2o[j] = row_index_mapping_n2o[j - 1];
            j = j - 1;
        }
        diag_entry[j] = lbuff0;
        row_index_mapping_n2o[j] = lbuff1;
    }
    // Old index to new one
    for(i = 0; i < _rows; i++)
        row_index_mapping_o2n[row_index_mapping_n2o[i]] = i;
    // Maximum number of columns in the sparse table
    max_columns = diag_entry[0];
    //--- End of sorting
    //
    //--- Create sparse table
    //
    num_column_entries = new long[max_columns];
    entry_column = new long[_size_entry_column];
    // 1. Count entries in each column in sparse table
    for (i = 0; i < max_columns; i++)
        num_column_entries[i] = 0;
    for (i = 0; i < _rows; i++)
        // 'diag_entry' still is used as a temporary array
        // it stores that numbers of nodes connect to this nodes
        for (j = 0; j < diag_entry[i]; j++)
            num_column_entries[j]++;
    // 2. Fill the sparse table, i.e. store all its entries to
    //    entry_column
    lbuff0 = 0;
    for (i = 0; i < max_columns; i++)
        for (j = 0; j < num_column_entries[i]; j++)
        {
            // ii is the real row index of this entry in matrix
            ii = row_index_mapping_n2o[j];
            // jj is the real column index of this entry in matrix
            jj = m_dom.get_node_conneted_nodes(ii, i);
            entry_column[lbuff0] = jj;
            // Till to this stage, 'diag_entry' is really used to store indices of the diagonal entries.
            // Hereby, 'index' refers to the index in entry_column array.
            if(ii == jj)
                diag_entry[ii] = lbuff0;
            //
            lbuff0++;
        }
}
/*\!
 ********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
void SparseTable::Write(std::ostream &os)
{
    long i, k, counter = 0;

    os.width(10);
    os << "Symmetry: " << symmetry << std::endl;
    os << "\n*** Row index  " << std::endl;

    if(storage_type == CRS)
    {
        os << "\n*** Sparse entry  " << std::endl;
        for (i = 0; i < _rows; i++)
        {
            for (k = num_column_entries[i]; k < num_column_entries[i + 1]; k++)
                os << entry_column[k] + 1 << " ";
            os << std::endl;
        }
    }
    else if(storage_type == JDS)
    {
        for (i = 0; i < _rows; i++)
            os << row_index_mapping_n2o[i] + 1 << std::endl;
        //
        os << "\n*** Sparse entry  " << std::endl;
        for (k = 0; k < max_columns; k++)
        {
            os << "--Column: " << k + 1 << std::endl;
            for (i = 0; i < num_column_entries[k]; i++)
            {
                os << entry_column[counter] + 1 << std::endl;
                counter++;
            }
            os << std::endl;
        }
    }
}

/*\!
 ********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
SparseTable::~SparseTable()
{
    if(entry_column)
        delete [] entry_column;
    if(num_column_entries)
        delete [] num_column_entries;
    if(row_index_mapping_n2o)
        delete [] row_index_mapping_n2o;
    if(row_index_mapping_o2n)
        delete [] row_index_mapping_o2n;
    if(diag_entry)
        delete [] diag_entry;
    entry_column = NULL;
    num_column_entries = NULL;
    row_index_mapping_n2o = NULL;
    row_index_mapping_o2n = NULL;
    diag_entry = NULL;
}
} //end
