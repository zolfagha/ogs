
#include <gtest/gtest.h>

#include <cmath>
#include <set>

#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/UnstructuredMesh.h"
#include "MeshLib/Tools/Mapping.h"
#include "MeshLib/Topology/Topology.h"
#include "MeshLib/Tools/MeshGenerator.h"

#include "TestUtil.h"

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

UnstructuredMesh * createLineElement(GeoLib::Point &pt1, GeoLib::Point &pt2, CoordinateSystemType::type coord_type)
{
    UnstructuredMesh *msh = new UnstructuredMesh(coord_type);
    msh->addNode(pt1);
    msh->addNode(pt2);
    Line *e = new Line();
    e->setNodeID(0,0);
    e->setNodeID(1,1);
    msh->addElement(e);

    return msh;
}

UnstructuredMesh * createTriangleElement(GeoLib::Point &pt1, GeoLib::Point &pt2, GeoLib::Point &pt3, CoordinateSystemType::type coord_type)
{
    UnstructuredMesh *msh = new UnstructuredMesh(coord_type);
    msh->addNode(pt1);
    msh->addNode(pt2);
    msh->addNode(pt3);
    Triangle *e = new Triangle();
    e->setNodeID(0,0);
    e->setNodeID(1,1);
    e->setNodeID(2,2);
    msh->addElement(e);

    return msh;
}

TEST(Mesh, MappingLineX1)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(1,0,0), GeoLib::Point(3,0,0), CoordinateSystemType::X);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineX2)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(3,0,0), GeoLib::Point(1,0,0), CoordinateSystemType::X);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineY1)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(0,1,0), GeoLib::Point(0,3,0), CoordinateSystemType::Y);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineY2)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(0,3,0), GeoLib::Point(0,1,0), CoordinateSystemType::Y);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineZ1)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(0,0,1), GeoLib::Point(0,0,3), CoordinateSystemType::Z);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineZ2)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(0,0,3), GeoLib::Point(0,0,1), CoordinateSystemType::Z);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineXY1)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(0,3,0), GeoLib::Point(0,1,0), CoordinateSystemType::XY);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineXY2)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(1,1,0), GeoLib::Point(2,2,0), CoordinateSystemType::XY);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(sqrt(2.),0,0), *p2);
}

TEST(Mesh, MappingLineXZ1)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(0,0,3), GeoLib::Point(0,0,1), CoordinateSystemType::XZ);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
}

TEST(Mesh, MappingLineXZ2)
{
    UnstructuredMesh* msh = createLineElement(GeoLib::Point(1,0,1), GeoLib::Point(2,0,2), CoordinateSystemType::XZ);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(sqrt(2.),0,0), *p2);
}

TEST(Mesh, MappingTriXY1)
{
    UnstructuredMesh* msh = createTriangleElement(GeoLib::Point(1,1,0), GeoLib::Point(3,1,0), GeoLib::Point(1,3,0), CoordinateSystemType::XY);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);
    GeoLib::Point* p3 = emap.getNodePoint(2);

    ASSERT_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_EQ(GeoLib::Point(2,0,0), *p2);
    ASSERT_EQ(GeoLib::Point(0,2,0), *p3);
}

TEST(Mesh, MappingTriXY2)
{
    UnstructuredMesh* msh = createTriangleElement(GeoLib::Point(3,1,0), GeoLib::Point(1,3,0), GeoLib::Point(1,1,0), CoordinateSystemType::XY);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);
    GeoLib::Point* p3 = emap.getNodePoint(2);

    ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point(sqrt(8.),0,0), *p2);
    ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point(sqrt(8.)*0.5,sqrt(2.),0), *p3);
}

TEST(Mesh, MappingTriXY3)
{
    UnstructuredMesh* msh = createTriangleElement(GeoLib::Point(1,3,0), GeoLib::Point(1,1,0), GeoLib::Point(3,1,0), CoordinateSystemType::XY);
    EleMapLocalCoordinates emap(msh, msh->getElemenet(0), msh->getGeometricProperty()->getCoordinateSystem());

    GeoLib::Point* p1 = emap.getNodePoint(0);
    GeoLib::Point* p2 = emap.getNodePoint(1);
    GeoLib::Point* p3 = emap.getNodePoint(2);

    ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point(0,0,0), *p1);
    ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point(2,0,0), *p2);
    ASSERT_DOUBLE_ARRAY_EQ(GeoLib::Point(2,2,0), *p3);
}
