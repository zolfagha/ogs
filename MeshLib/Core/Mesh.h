
#pragma once

#include <cmath>
#include <vector>
#include "MemoryTools.h"
#include "Node.h"
#include "Element.h"
#include "ElementFactory.h"

//------------------------------------------------------------------------
namespace MeshLib
{
////-----------------------------------------------------------------------------
//// Mesh
////-----------------------------------------------------------------------------
//class IMesh
//{
//public:
//    virtual ~IMesh(){};
//
//    virtual size_t getNumberOfNodes() const = 0;
//    virtual size_t getNumberOfElements() const = 0;
//    virtual void getNodeCoordinates( size_t node_id, double pt[3] ) const = 0;
//    virtual IElement* getElemenet( size_t element_id ) const = 0;
//
//    virtual Node* getNode( size_t id ) const = 0;
//};
//
//class StructuredMesh : public IMesh
//{
//private:
//    double _origin[3];
//    double _length[3];
//    double  _unit_length[3];
//    size_t  _number_of_nodes_per_dimension[3];
//public:
//};
//
//class UnstructuredMesh : public IMesh
//{
//private:
//    std::vector<Node*> _list_nodes;
//    std::vector<IElement*> _list_elements;
//
//public:
//    UnstructuredMesh(){};
//    virtual ~UnstructuredMesh(){
//        destroyStdVectorWithPointers(_list_nodes);
//        destroyStdVectorWithPointers(_list_elements);
//    };
//
//    virtual size_t getNumberOfNodes() const { return _list_nodes.size(); };
//    virtual size_t getNumberOfElements() const { return _list_elements.size(); };
//    virtual void getNodeCoordinates(size_t node_id, double pt[3]) const {
//        Node* nod = _list_nodes[node_id];
//        pt[0] = (*nod)[0];
//        pt[1] = (*nod)[1];
//        pt[2] = (*nod)[2];
//    }
//
//    size_t setNode( size_t node_id, double x, double y, double z ) {
//        size_t new_node_id = node_id;
//        if (node_id<_list_nodes.size()) {
//            Node *node = this->_list_nodes.at(node_id);
//            (*node)[0] = x;
//            (*node)[1] = y;
//            (*node)[2] = z;
//        } else {
//            new_node_id = this->_list_nodes.size();
//            this->_list_nodes.push_back(new Node(new_node_id, x, y, z));
//        }
//        return new_node_id;
//    };
//
//    size_t addElement( IElement *e) {
//        e->setElementID(_list_elements.size());
//        _list_elements.push_back(e);
//        return e->getElementID();
//    };
//
//    virtual IElement* getElemenet( size_t element_id ) const {
//        assert(element_id<_list_elements.size());
//        return _list_elements.at(element_id);
//    };
//
//    virtual Node* getNode( size_t id ) const {
//        return _list_nodes.at(id);
//    };
//
//    void construct() {
//
//        //set node connectivity
//        const size_t n_ele = this->getNumberOfElements();
//        for (size_t i=0; i<n_ele; i++) {
//            const IElement* e = _list_elements[i];
//            for (size_t j=0; j<e->getNumberOfNodes(); j++) {
//                Node* nod_j = getNode(e->getNodeID(j));
//                for (size_t k=j+1; k<e->getNumberOfNodes(); k++) {
//                    Node* nod_k = getNode(e->getNodeID(k));
//                    nod_j->addConnectedNode(e->getNodeID(k));
//                    nod_k->addConnectedNode(e->getNodeID(j));
//                }
//            };
//        }
//    };
//};

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

