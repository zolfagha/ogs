#pragma  once

#include <vector>
#include "Base/MemoryTools.h"
#include "IElement.h"
#include "Node.h"

namespace MeshLib
{

/**
 * 
 */
template <ElementType::type TYPE, size_t NUMBER_OF_NODES, size_t DIMENSION, size_t NUMBER_OF_FACES, size_t NUMER_OF_EDGES>
class TemplateUnstructuredElement : public IElement
{
private:
    INode* _list_nodes[NUMBER_OF_NODES];
    size_t _list_node_id[NUMBER_OF_NODES];
    IElementMapping *_geo_map;
    std::vector<IElement*> _list_edges;

    IElement* createEdgeElement(size_t edge_id);
public:
    TemplateUnstructuredElement() {
        initialize();
    };
    TemplateUnstructuredElement(size_t element_id) {
        initialize();
        this->setElementID(element_id);
    }
    virtual ~TemplateUnstructuredElement() 
    {
        destroyStdVectorWithPointers(_list_edges);
    };

    virtual void initialize() {
        for (int i=0; i<NUMBER_OF_NODES; i++)
            _list_node_id[i] = 0;
        for (int i=0; i<NUMER_OF_EDGES; i++)
            _list_edges[i] = 0;
        _geo_map = 0;
        if (NUMER_OF_EDGES>0)
            _list_edges.resize(NUMER_OF_EDGES, 0);
    }

    virtual ElementType::type getElementType() const {return TYPE;};
    virtual size_t getNumberOfNodes() const {return NUMBER_OF_NODES;};
    virtual size_t getDimension() const {return DIMENSION;};
    virtual size_t getNumberOfFaces() const {return NUMBER_OF_FACES;};
    virtual size_t getNumberOfEdges() const {return NUMER_OF_EDGES;};
    virtual IElement* getEdgeElement(size_t edge_id) 
    {
        if (_list_edges[edge_id]==0)
            _list_edges[edge_id] = createEdgeElement(edge_id);
        return _list_edges[edge_id];
    }

    virtual void setNodeID(size_t local_node_id, size_t node_id) {
        assert (local_node_id < NUMBER_OF_NODES);
        _list_node_id[local_node_id] = node_id;
    }
    virtual size_t getNodeID(size_t local_node_id) const {
        assert (local_node_id < NUMBER_OF_NODES);
        return _list_node_id[local_node_id];
    };
    void setNode(size_t local_node_id, INode* nod) {
        _list_nodes[local_node_id] = nod;
    }
    INode* getNode(size_t local_node_id) const {
        return _list_nodes[local_node_id];
    }

    virtual const GeoLib::Point* getNodeLocalCoordinates( size_t i_nod ) const 
    {
        return _list_nodes[i_nod]->getData();
    };

    virtual IElementMapping* getMappedGeometry() {
        if (_geo_map==0)
            _geo_map = new EleMapLocalCoordinates(this);
        return _geo_map;
    };
};

// elements
typedef TemplateUnstructuredElement<ElementType::LINE, 2, 1, 2, 0> Line;
typedef TemplateUnstructuredElement<ElementType::TRIANGLE, 3, 2, 3, 3> Triangle;
typedef TemplateUnstructuredElement<ElementType::QUAD, 4, 2, 4, 4> Quadrirateral;
typedef TemplateUnstructuredElement<ElementType::TETRAHEDRON, 4, 3, 4, 6> Tetrahedron;
typedef TemplateUnstructuredElement<ElementType::PYRAMID, 5, 3, 5, 8> Pyramid;
typedef TemplateUnstructuredElement<ElementType::PRISM, 6, 3, 5, 9> Prism;
typedef TemplateUnstructuredElement<ElementType::HEXAHEDRON, 8, 3, 6, 12> Hexahedron;

template <> 
IElement* Line::createEdgeElement(size_t edge_id)
{
    return 0;
}

template <> 
IElement* Triangle::createEdgeElement(size_t edge_id)
{
    Line *e = new Line();
    e->setNodeID(0, e->getNodeID(edge_id));
    e->setNodeID(1, e->getNodeID((edge_id+1)%3));
    return e;
}

template <> 
IElement* Quadrirateral::createEdgeElement(size_t edge_id)
{
    Line *e = new Line();
    e->setNodeID(0, e->getNodeID(edge_id));
    e->setNodeID(1, e->getNodeID((edge_id+1)%4));
    return e;
}

} // end namespace

