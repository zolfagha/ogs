
#pragma once

#include <set>

#include "MathLib/LinAlg/Sparse/SparseTableCRS.h"

#include "MeshLib/Topology/Topology.h"

#include "DomainDecomposition.h"
#include "DofMapManager.h"

namespace DiscreteLib
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

    //static void createRowMajorSparsityFromNodeConnectivity(const MeshLib::ITopologyNode2Nodes &topo_node2nodes, const DofMap &dof, MathLib::RowMajorSparsity &row_major_entries)
    //{
    //    const size_t n_nodes = topo_node2nodes.getNumberOfNodes();
    //    row_major_entries.resize(dof.getNumberOfActiveDoFs());

    //    size_t i_rows = 0;
    //    for (size_t i=0; i<n_nodes; i++) {
    //        if (!dof.isActiveDoF(i)) continue;
    //        // search connected nodes
    //        std::set<size_t> &setConnection = row_major_entries[i_rows++];
    //        setConnection.insert(dof.getEqsID(i));
    //        const std::set<size_t> connected_nodes = topo_node2nodes.getConnectedNodes(i);
    //        for (std::set<size_t>::const_iterator it=connected_nodes.begin(); it!=connected_nodes.end(); it++) {
    //            if (dof.isActiveDoF(*it)) setConnection.insert(dof.getEqsID(*it));
    //        }
    //    }
    //}

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

class SparsityBuilderDummy
{
public:
    SparsityBuilderDummy(MeshLib::IMesh&, DofEquationIdTable&, MathLib::RowMajorSparsity&) {};
};

class SparsityBuilderFromNodeConnectivity
{
public:
    SparsityBuilderFromNodeConnectivity(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
    {
        MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(&msh);
        if (dofManager.getNumberOfVariables()==1) {
            SparsityBuilder::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparse);
        } else {
            SparsityBuilder::createRowMajorSparsityForMultipleDOFs(topo_node2nodes, dofManager.getNumberOfVariables(), sparse);
        }
    }
};

class SparsityBuilderFromNodeConnectivityWithGhostDoFs
{
public:
    SparsityBuilderFromNodeConnectivityWithGhostDoFs(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
    {
        MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(&msh);
        if (dofManager.getNumberOfVariables()==1) {
            SparsityBuilder::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparse);
        } else {
            SparsityBuilder::createRowMajorSparsityForMultipleDOFs(topo_node2nodes, dofManager.getNumberOfVariables(), sparse);
        }
        size_t n_rows = dofManager.getTotalNumberOfActiveDoFsWithoutGhost();
        sparse.erase(sparse.begin()+n_rows, sparse.end());
    }
};

//class SparsityBuilderFromNodeConnectivityWithInactiveDoFs
//{
//public:
//    SparsityBuilderFromNodeConnectivityWithInactiveDoFs(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
//    {
//        MeshLib::TopologyNode2NodesConnectedByElements topo_node2nodes(&msh);
//        if (dofManager.getNumberOfVariables()==1) {
//            SparsityBuilder::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, *dofManager.getVariableDoF(0), sparse);
//        } else {
//            SparsityBuilder::createRowMajorSparsityForMultipleDOFs(topo_node2nodes, dofManager.getNumberOfVariables(), sparse);
//        }
//    }
//};

class SparsityBuilderFromLocalSparsity
{
public:
    SparsityBuilderFromLocalSparsity(std::vector<MathLib::RowMajorSparsity*> &list_local_sparse, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
    {

    }
};

template <class T_SPARSITY_BUILDER>
class SparsityBuilderFromDDC
{
public:
    SparsityBuilderFromDDC(DDCGlobal &ddc_global, DofEquationIdTable &dofManager, MathLib::RowMajorSparsity &sparse)
    {
        const size_t n_dom = ddc_global.getNumberOfSubDomains();
        const size_t n_dofs = dofManager.getTotalNumberOfActiveDoFs();
        sparse.resize(n_dofs);

        std::vector<MathLib::RowMajorSparsity> local_sparse(ddc_global.getNumberOfSubDomains());
        for (size_t i=0; i<n_dom; i++) {
            DDCSubDomain* dom = ddc_global.getSubDomain(i);
            T_SPARSITY_BUILDER builder(*dom->getLoalMesh(), dofManager, local_sparse[i]);
        }

        for (size_t i=0; i<n_dom; i++) {
            DDCSubDomain* dom = ddc_global.getSubDomain(i);
            IDDCGlobaLocalMapping* mapping = dom->getGlobalLocalIdMap();
            MathLib::RowMajorSparsity &local_sp = local_sparse[i];
            for (size_t j=0; j<local_sp.size(); j++) {
                const std::set<size_t> &local_row = local_sp[j];
                const size_t i_gloal = mapping->local2global(j);
                std::set<size_t> &global_row = sparse[i_gloal];
                for (std::set<size_t>::iterator itr=local_row.begin(); itr!=local_row.end(); ++itr) {
                    global_row.insert(mapping->local2global(*itr));
                }
            }
        }

        //size_t i_row = 0;
        //for (size_t i=0; i<local_sparse.size(); i++) {
        //    size_t offset = 0;
        //    MathLib::RowMajorSparsity &local_sp = local_sparse[i];
        //    for (size_t j=0; j<local_sp.size(); j++) {
        //        std::set<size_t> &setConnection = sparse[i_row++];
        //        std::set<size_t> &local_conn = local_sp[j];
        //        for (std::set<size_t>::iterator itr = local_conn.begin(); itr!=local_conn.end(); ++itr) {
        //            setConnection.insert(offset + *itr);
        //        }
        //    }
        //    offset += local_sp.size();
        //}
    }
};

}
