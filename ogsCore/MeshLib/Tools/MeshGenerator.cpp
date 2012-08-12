/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshGenerator.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "MeshGenerator.h"

#include <memory>
#include "GeoLib/Point.h"
#include "MeshLib/Topology/Topology.h"

namespace MeshLib
{

MeshLib::UnstructuredMesh* MeshGenerator::generateLineMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
{
    MeshLib::UnstructuredMesh* msh = new MeshLib::UnstructuredMesh(MeshLib::CoordinateSystemType::X);

    const double unit_length = length / subdivision;
    const size_t n_nodes_per_axis = subdivision+1;
    const size_t n_eles = subdivision;
    //const size_t n_nodes = subdivision+1;

    //nodes
    size_t node_id(0);
    for (size_t i_z=0; i_z<n_nodes_per_axis; i_z++) {
        const double x = unit_length*i_z;
        GeoLib::Point p(x+origin_x, origin_y, origin_z);
        msh->setNodeCoordinates(node_id++, p);
        //msh->setNodeCoordinates(node_id++, GeoLib::Point(x+origin_x, origin_y, origin_z));
    }

    //elements
    //size_t ele_id(0);
    for (size_t i_z=0; i_z<n_eles; i_z++) {
        Line *e = new Line();
        e->setNodeID(0, i_z);
        e->setNodeID(1, i_z+1);
        msh->addElement(e);
    }

    return msh;
}

MeshLib::UnstructuredMesh* MeshGenerator::generateRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
{

    MeshLib::UnstructuredMesh* msh = new MeshLib::UnstructuredMesh(MeshLib::CoordinateSystemType::XY);

    //size_t n_eles = static_cast<size_t>(pow(static_cast<double>(subdivision), 2));
    //size_t n_nodes = static_cast<size_t>(pow(static_cast<double>(subdivision+1), 2));
    const double unit_length = length / subdivision;
    const size_t n_nodes_per_axis = subdivision+1;

    //nodes
    size_t node_id(0);
    const double z = origin_z;
    for (size_t j_y=0; j_y<n_nodes_per_axis; j_y++) {
        const double y = unit_length*j_y + origin_y;
        for (size_t k_x=0; k_x<n_nodes_per_axis; k_x++) {
            const double x = unit_length*k_x + origin_x;
            GeoLib::Point p(x, y, z);
            msh->setNodeCoordinates(node_id++, p);
            //msh->setNodeCoordinates(node_id++, GeoLib::Point(x, y, z));
        }
    }

    //elements
    //size_t ele_id(0);
    for (size_t j=0; j<subdivision; j++) {
        const size_t offset_y1 = j*n_nodes_per_axis;
        const size_t offset_y2 = (j+1)*n_nodes_per_axis;
        for (size_t k=0; k<subdivision; k++) {
            Quadrirateral *e = new Quadrirateral();
            e->setNodeID(0, offset_y1+k);
            e->setNodeID(1, offset_y1+k+1);
            e->setNodeID(2, offset_y2+k+1);
            e->setNodeID(3, offset_y2+k);
            //for (size_t l=0; l<4; l++)
            //    e->setNode(l, msh->getNode(e->getNodeID(l)));
            msh->addElement(e);
        }
    }

    return msh;

};

StructuredMesh<ElementShape::QUAD>* MeshGenerator::generateStructuredRegularQuadMesh(const double length, const size_t subdivision, const double origin_x, const double origin_y, const double origin_z) 
{
    GeoLib::Point org(origin_x, origin_y, origin_z);
    GeoLib::Point pt_length(length, length, .0);
    size_t unit_div[3] = {subdivision, subdivision, 1};

    return new StructuredMesh<ElementShape::QUAD>(MeshLib::CoordinateSystemType::XY, org, pt_length.getData(), unit_div);

};

void MeshGenerator::generateSubMesh(const MeshLib::IMesh &src, const std::vector<size_t> &list_e, MeshLib::IMesh* &dest, BaseLib::BidirectionalMap<size_t, size_t> &map_global2local)
{
    MeshLib::UnstructuredMesh* new_msh = new MeshLib::UnstructuredMesh(src.getGeometricProperty()->getCoordinateSystem()->getType());
    std::set<size_t> list_nodes_subset;
    for (size_t i=0; i<list_e.size(); i++) {
        const MeshLib::IElement *e = src.getElemenet(list_e[i]);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            list_nodes_subset.insert(e->getNodeID(j));
        }
    }

    for (std::set<size_t>::iterator itr=list_nodes_subset.begin(); itr!=list_nodes_subset.end(); ++itr) {
        GeoLib::Point p(src.getNodeCoordinates(*itr));
        size_t new_node_id = new_msh->addNode(p);
        //size_t new_node_id = new_msh->addNode(src.getNodeCoordinates(*itr));
        map_global2local.insert(*itr, new_node_id);
    }

    for (size_t i=0; i<list_e.size(); i++) {
        const MeshLib::IElement *e = src.getElemenet(list_e[i]);
        MeshLib::IElement *new_e = e->clone();
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            new_e->setNodeID(j, map_global2local.mapAtoB(e->getNodeID(j)));
        }
        new_msh->addElement(new_e);
    }
    dest = new_msh;
}

}
