/**
 * \file msh_mesh.h
 */

#ifndef msh_mesh_INC
#define msh_mesh_INC

#include <cmath>
#include <vector>
#include "MemoryTools.h"
#include "Point.h"

//------------------------------------------------------------------------
namespace MeshLib
{

//-----------------------------------------------------------------------------
// Node
//-----------------------------------------------------------------------------
class Node : public GEOLIB::Point
{
private:
    size_t _node_id;
public:
    Node (size_t id, double x, double y, double z) {
        this->_node_id = id;
        this->_x[0] = x;
        this->_x[1] = y;
        this->_x[2] = z;
    };


};

//-----------------------------------------------------------------------------
// Element
//-----------------------------------------------------------------------------

struct ElementType
{
    enum type {
        LINE = 1,
        QUAD = 2,
        HEXAHEDRON = 3,
        TRIANGLE = 4,
        TETRAHEDRON = 5,
        PRISM = 6,
        PYRAMID = 7,
        INVALID = -1
    };
};

class IElement
{
private:
    size_t _element_id;
    size_t _group_id;
public:
    IElement():_element_id(0), _group_id(0) {};
    virtual ~IElement() {};

    void setElementID(size_t id) {_element_id = id;};
    size_t getElementID() const {return _element_id;};
    size_t setGroupID(size_t id) {_group_id = id;};
    size_t getGroupID() const {return _group_id;};

    virtual ElementType::type getElementType() const = 0;
    virtual size_t getElementDimension() const = 0;
    virtual size_t getNumberOfNodes() const = 0;
    virtual size_t getNumberOfFaces() const = 0;
    virtual size_t getNumberOfEdges() const = 0;

    virtual void setNodeID(size_t local_node_id, size_t node_id) = 0;
    virtual size_t getNodeID(size_t local_node_id) const = 0;
};

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

//-----------------------------------------------------------------------------
// Mesh
//-----------------------------------------------------------------------------
class IMesh
{
public:
    virtual ~IMesh(){};

    virtual size_t getNumberOfNodes() const = 0;
    virtual size_t getNumberOfElements() const = 0;
    virtual void getNodeCoordinates( size_t node_id, double pt[3] ) const = 0;
    virtual IElement* getElemenet( size_t element_id ) const = 0;
};

class StructuredMesh : public IMesh
{
private:
    double _origin[3];
    double _length[3];
    double  _unit_length[3];
    size_t  _number_of_nodes_per_dimension[3];
public:
};

class UnstructuredMesh : public IMesh
{
private:
    std::vector<Node*> _list_nodes;
    std::vector<IElement*> _list_elements;

public:
    UnstructuredMesh(){};
    virtual ~UnstructuredMesh(){
        destroyStdVectorWithPointers(_list_nodes);
        destroyStdVectorWithPointers(_list_elements);
    };

    virtual size_t getNumberOfNodes() const { return _list_nodes.size(); };
    virtual size_t getNumberOfElements() const { return _list_elements.size(); };
    virtual void getNodeCoordinates(size_t node_id, double pt[3]) const {
        Node* nod = _list_nodes[node_id];
        pt[0] = (*nod)[0];
        pt[1] = (*nod)[1];
        pt[2] = (*nod)[2];
    }

    size_t setNode( size_t node_id, double x, double y, double z ) {
        size_t new_node_id = node_id;
        if (node_id<_list_nodes.size()) {
            Node *node = this->_list_nodes.at(node_id);
            (*node)[0] = x;
            (*node)[1] = y;
            (*node)[2] = z;
        } else {
            new_node_id = this->_list_nodes.size();
            this->_list_nodes.push_back(new Node(new_node_id, x, y, z));
        }
        return new_node_id;
    };

    size_t addElement( IElement *e) {
        e->setElementID(_list_elements.size());
        _list_elements.push_back(e);
        return e->getElementID();
    };

    virtual IElement* getElemenet( size_t element_id ) const {
        assert(element_id<_list_elements.size());
        return _list_elements.at(element_id);
    }
};

template <size_t N_DIM>
class StaticStructuredMesh : public IMesh
{
private:
    size_t  _dimensions[N_DIM];
    double  _origin[N_DIM];
    size_t  _number_of_nodes[N_DIM];
    double  _unit_length[N_DIM];
public:
};

template <size_t  N_DIM, size_t  N_NODES, size_t  N_ELEMENTS>
class StaticUnstructuredMesh : public IMesh
{
private:
    size_t  _nodes[N_NODES];
    size_t  _elements[N_ELEMENTS];
    size_t  _edges;
    size_t  _faces;

public:
    StaticUnstructuredMesh(){};
    virtual ~StaticUnstructuredMesh(){};

    virtual size_t getNumberOfNodes() const { return N_NODES; };
    virtual size_t getNumberOfElements() const { return N_ELEMENTS; };
};

class HierarchicalMesh : public IMesh
{
private:
};

} // end namespace

#endif
