
#include <gtest/gtest.h>

#include <set>

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Topology/Topology.h"
#include "MeshLib/Tools/MeshGenerator.h"

using namespace MeshLib;

TEST(Mesh, topoN2N)
{
    IMesh* msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    TopologyNode2Nodes topo_node2nodes(msh);

    ASSERT_EQ(topo_node2nodes.getNumberOfNodes(), msh->getNumberOfNodes());
    size_t expected_conn_nodes0[] = {1,3};
    size_t expected_conn_nodes1[] = {0,2,4};
    size_t expected_conn_nodes4[] = {1,3,5,7};
    ASSERT_EQ(std::set<size_t>(expected_conn_nodes0, expected_conn_nodes0+2), topo_node2nodes.getConnectedNodes(0));
    ASSERT_EQ(std::set<size_t>(expected_conn_nodes1, expected_conn_nodes1+3), topo_node2nodes.getConnectedNodes(1));
    ASSERT_EQ(std::set<size_t>(expected_conn_nodes4, expected_conn_nodes4+4), topo_node2nodes.getConnectedNodes(4));
}
