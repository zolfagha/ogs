
#pragma once

#include <vector>
#include "IMesh.h"
#include "IElement.h"
#include "Node.h"
#include "Base/MemoryTools.h"

namespace MeshLib
{

template<typename Tpos, size_t N_DIM>
class UnstructuredMesh : public IMesh<Tpos>
{
private:
    std::vector<Node<Tpos>*> _list_nodes;
    std::vector<IElement*> _list_elements;

public:
    UnstructuredMesh()
    {
    };
    virtual ~UnstructuredMesh(){
        destroyStdVectorWithPointers(_list_nodes);
        destroyStdVectorWithPointers(_list_elements);
    };

    virtual size_t getDimension() const
    {
      return N_DIM;
    }

    virtual size_t getNumberOfNodes() const { return _list_nodes.size(); };
    virtual size_t getNumberOfElements() const { return _list_elements.size(); };
    virtual void getNodeCoordinates(size_t node_id, Tpos &pt) const {
        Node<Tpos>* nod = _list_nodes[node_id];
        pt = *nod->getData();
    }

    size_t setNode( size_t node_id, Tpos &x ) {
        size_t new_node_id = node_id;
        if (node_id<_list_nodes.size()) {
            Node<Tpos> *node = this->_list_nodes.at(node_id);
            node->setX(x);
        } else {
            new_node_id = this->_list_nodes.size();
            this->_list_nodes.push_back(new Node<Tpos>(new_node_id, x));
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
    };

    virtual INode<Tpos>* getNode( size_t id ) const {
        return _list_nodes.at(id);
    };

    void construct() {

        //set node connectivity
        const size_t n_ele = this->getNumberOfElements();
        for (size_t i=0; i<n_ele; i++) {
            const IElement* e = _list_elements[i];
            for (size_t j=0; j<e->getNumberOfNodes(); j++) {
                Node* nod_j = getNode(e->getNodeID(j));
                for (size_t k=j+1; k<e->getNumberOfNodes(); k++) {
                    Node* nod_k = getNode(e->getNodeID(k));
                    nod_j->addConnectedNode(e->getNodeID(k));
                    nod_k->addConnectedNode(e->getNodeID(j));
                }
            };
        }
    };
};

}
