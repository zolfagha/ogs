
#pragma once

#include "IElement.h"

namespace MeshLib
{

template <ElementShape::type TYPE, size_t DIMENSION, size_t NUMBER_OF_FACES, size_t NUMER_OF_EDGES, size_t NUMBER_OF_NODES1, size_t NUMBER_OF_NODES2>
class TemplateElementTopology
{
public:
    static ElementShape::type getShapeType() {return TYPE;};
    static size_t getDimension() {return DIMENSION;};
    static size_t getNumberOfFaces() {return NUMBER_OF_FACES;};
    static size_t getNumberOfEdges() {return NUMER_OF_EDGES;};
    static size_t getNumberOfNodes(size_t order)
    {
        switch (order) {
        case 1: return NUMBER_OF_NODES1;
        case 2: return NUMBER_OF_NODES2;
        }
        return 0;
    }
};

class LineTopology : public TemplateElementTopology<ElementShape::LINE, 1, 0, 1, 2, 3>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class TriangleTopology : public TemplateElementTopology<ElementShape::TRIANGLE, 2, 0, 3, 3, 6>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class Quad8Topology : public TemplateElementTopology<ElementShape::QUAD, 2, 0, 4, 4, 8>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class Quad9Topology : public TemplateElementTopology<ElementShape::QUAD, 2, 0, 4, 4, 9>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
};

class TetraTopology : public TemplateElementTopology<ElementShape::TETRAHEDRON, 3, 4, 6, 4, 10>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t face_id) { return ElementShape::TRIANGLE; }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

class HexTopology : public TemplateElementTopology<ElementShape::HEXAHEDRON, 3, 6, 12, 8, 20>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t face_id) { return ElementShape::QUAD; }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

class PrismTopology : public TemplateElementTopology<ElementShape::PRISM, 3, 5, 9, 6, 15>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t face_id)
    { 
        if (face_id<2) return ElementShape::TRIANGLE;
        else return ElementShape::QUAD;
    }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

class PyramidTopology : public TemplateElementTopology<ElementShape::PYRAMID, 3, 5, 8, 5, 13>
{
public:
    static ElementShape::type getEdgeElementType(size_t edge_id) { return ElementShape::LINE; }
    static void getLocalNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_local_node_ids);
    static ElementShape::type getFaceElementType(size_t face_id)
    { 
        if (face_id==0) return ElementShape::QUAD;
        else return ElementShape::TRIANGLE;
    }
    static void getLocalNodeIDsOfFaceElement(size_t order, size_t face_id, std::vector<size_t> &vec_local_node_ids);
};

} //end namespace

