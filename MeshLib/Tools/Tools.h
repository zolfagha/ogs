
#pragma once

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Tools/MeshNodesAlongPolyline.h"

#include <vector>
#include <algorithm>

namespace MeshLib
{

/// seeks all nodes located on a given polyline and returns a list of the found node pointers.
void findNodesOnPolyline(IMesh const* msh, GeoLib::Polyline const* poly, std::vector<INode*> *vec_nodes)
{
    MeshNodesAlongPolyline obj(poly, (MeshLib::UnstructuredMesh2d*)msh);
    std::vector<size_t> vec_node_id = obj.getNodeIDs();
    vec_nodes->resize(vec_node_id.size());
    for (size_t i=0; i<vec_node_id.size(); i++) {
        (*vec_nodes)[i] = msh->getNode(vec_node_id[i]);
    }
};

///
void findNodesOnGeometry(IMesh const* msh, GeoLib::GeoObject const* obj, std::vector<INode*> *vec_nodes)
{
    switch (obj->getGeoType()) {
        case GeoLib::GeoObjType::POLYLINE:
            findNodesOnPolyline(msh, static_cast<GeoLib::Polyline const*>(obj), vec_nodes);
            break;
        default:
            throw std::exception("This geo type is not supported in MeshLib::findNodesOnGeometry()");
            break;
    }
};

/// get a list of elements connected to one of give nodes
void findConnectedElements(IMesh const* msh, const std::vector<INode*> &nodes, std::vector<IElement*> &connected_elements)
{
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        IElement* e = msh->getElemenet(i);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (std::find(nodes.begin(), nodes.end(), msh->getNode(e->getNodeID(j)))!=nodes.end()) {
                connected_elements.push_back(e);
                break;
            }
        }
    }
};

/// get a list of edge elements of given elements
void createEdgeElements(IMesh const* msh, const std::vector<IElement*> &selected_ele, std::vector<IElement*> &edges)
{
    for (size_t i=0; i<selected_ele.size(); i++) {
        IElement *e = selected_ele[i];
        for (size_t j=0; j<e->getNumberOfEdges(); j++) {
            IElement *edge = e->getEdgeElement(j);
            edges.push_back(edge);
        }
    }
};

/// find edge elements on a given polyline
void findEdgeElementsOnPolyline(IMesh const* msh, GeoLib::Polyline const* poly, std::vector<IElement*> *vec_edges_on_poly)
{
    // get a list of nodes on the polyline
    std::vector<INode*> nodes_on_poly;
    findNodesOnPolyline(msh, poly, &nodes_on_poly);
    // get a list of elements having the nodes
    std::vector<IElement*> elements_near_poly;
    findConnectedElements(msh, nodes_on_poly, elements_near_poly);
    // get a list of edges made of the nodes
    std::vector<IElement*> selected_edges;
    createEdgeElements(msh, elements_near_poly, selected_edges);
    const size_t nr_edges = selected_edges.size();
    for (size_t i=0; i<nr_edges; i++) {
        IElement *edge_e = selected_edges[i];
        // check
        size_t cnt_match = 0;
        for (size_t j=0; j<edge_e->getNumberOfNodes(); j++) {
            const INode* edge_nod = msh->getNode(edge_e->getNodeID(j));
            if (std::find(nodes_on_poly.begin(), nodes_on_poly.end(), edge_nod) == nodes_on_poly.end())
                cnt_match++;
        }
        // update the list
        if (cnt_match==edge_e->getNumberOfNodes())
            vec_edges_on_poly->push_back(edge_e);
    }
};

/// find boundary elements located on given geometries
void findBoundaryElementsOnGeometry(IMesh const* msh, GeoLib::GeoObject const* obj, std::vector<IElement*> *vec_eles)
{
    switch (obj->getGeoType()) {
        case GeoLib::GeoObjType::POLYLINE:
            findEdgeElementsOnPolyline(msh, static_cast<GeoLib::Polyline const*>(obj), vec_eles);
            break;
        default:
            throw std::exception("This geo type is not supported in MeshLib::findNodesOnGeometry()");
            break;
    }
};

}