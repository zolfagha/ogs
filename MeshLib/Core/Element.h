#pragma  once


#include "Node.h"

//------------------------------------------------------------------------
namespace MeshLib
{

//-----------------------------------------------------------------------------
// Element
//-----------------------------------------------------------------------------


template <ElementType::type TYPE, size_t NUMBER_OF_NODES, size_t DIMENSION, size_t NUMBER_OF_FACES, size_t NUMER_OF_EDGES>
class TemplateElement : public IElement
{
private:
    size_t _list_node_id[NUMBER_OF_NODES];
public:
    TemplateElement() {
        initialize();
    };
    TemplateElement(size_t element_id) {
        initialize();
        this->setElementID(element_id);
    }
    virtual ~TemplateElement() {};

    virtual void initialize() {
        for (int i=0; i<NUMBER_OF_NODES; i++)
            _list_node_id[i] = 0;
    }

    virtual ElementType::type getElementType() const {return TYPE;};
    virtual size_t getNumberOfNodes() const {return NUMBER_OF_NODES;};
    virtual size_t getElementDimension() const {return DIMENSION;};
    virtual size_t getNumberOfFaces() const {return NUMBER_OF_FACES;};
    virtual size_t getNumberOfEdges() const {return NUMER_OF_EDGES;};

    virtual void setNodeID(size_t local_node_id, size_t node_id) {
        assert (local_node_id < NUMBER_OF_NODES);
        _list_node_id[local_node_id] = node_id;
    }
    virtual size_t getNodeID(size_t local_node_id) const {
        assert (local_node_id < NUMBER_OF_NODES);
        return _list_node_id[local_node_id];
    };
};

// elements
typedef TemplateElement<ElementType::LINE, 2, 1, 2, 0> Line;
typedef TemplateElement<ElementType::TRIANGLE, 3, 2, 3, 3> Triangle;
typedef TemplateElement<ElementType::QUAD, 4, 2, 4, 4> Quadrirateral;
typedef TemplateElement<ElementType::TETRAHEDRON, 4, 3, 4, 6> Tetrahedron;
typedef TemplateElement<ElementType::PYRAMID, 5, 3, 5, 8> Pyramid;
typedef TemplateElement<ElementType::PRISM, 6, 3, 5, 9> Prism;
typedef TemplateElement<ElementType::HEXAHEDRON, 8, 3, 6, 12> Hexahedron;

} // end namespace

