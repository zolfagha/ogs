
#pragma once

#include <set>

#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

#include "MeshLib/Topology/Topology.h"

#include "NumLib/Discrete/DoF.h"

namespace NumLib
{

class SparsityBuilder
{
public:
    /**
     * create row-major sparsity based on node connectivities
     */
    static void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, MathLib::RowMajorSparsity &row_major_entries)
    {
        const size_t n_nodes = topo_node2nodes.getNumberOfNodes();
        row_major_entries.resize(n_nodes);

        for (size_t i=0; i<n_nodes; i++) {
            // search connected nodes
            std::set<size_t> &setConnection = row_major_entries[i];
            setConnection.insert(i);
            const std::set<size_t> connected_nodes = topo_node2nodes.getConnectedNodes(i);
            for (std::set<size_t>::const_iterator it=connected_nodes.begin(); it!=connected_nodes.end(); it++) {
                setConnection.insert(*it);
            }
        }
    }

    /**
     * create row-major sparsity for multiple DOFs with a single mesh
     */
    static void createRowMajorSparsityForMultipleDOFs(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, size_t n_dof, MathLib::RowMajorSparsity &row_major_entries)
    {
        const size_t n_nodes = topo_node2nodes.getNumberOfNodes();
        row_major_entries.resize(n_nodes*n_dof);

        MathLib::RowMajorSparsity sparsity_msh;
        createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparsity_msh);

        for (size_t i_dof=0; i_dof<n_dof; i_dof++) {
            for (size_t i=0; i<n_nodes; i++) {
                std::set<size_t> &msh_row_entries = sparsity_msh[i];
                std::set<size_t> &row_entries = row_major_entries[i];
                for (std::set<size_t>::iterator it=msh_row_entries.begin(); it!=msh_row_entries.end(); ++it) {
                    size_t nod_id = *it;
                    for (size_t j=0; j<n_dof; j++) {
                        size_t eqs_id = 0;
                        if (0) {
                            eqs_id = nod_id*n_dof + j;
                        } else {
                            eqs_id = i + j*n_nodes;
                        }
                        row_entries.insert(eqs_id);
                    }
                }
            }
        }
    }
};



}
