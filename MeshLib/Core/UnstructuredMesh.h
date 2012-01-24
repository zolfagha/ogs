
#pragma once

#include <iostream>
#include <vector>

#include "Base/MemoryTools.h"
#include "MathLib/Vector.h"

#include "IMesh.h"
#include "IElement.h"
#include "INode.h"
#include "Node.h"

namespace MeshLib
{

template<size_t N_DIM, CoordinateSystem::CoordinateSystemType coord_type>
class TemplateUnstructuredMesh : public IMesh
{
private:
    std::vector<Node*> _list_nodes;
    std::vector<IElement*> _list_elements;
    CoordinateSystem _coord_system;

public:
    TemplateUnstructuredMesh() : _coord_system(coord_type)
    {
    };
    virtual ~TemplateUnstructuredMesh(){
        destroyStdVectorWithPointers(_list_nodes);
        destroyStdVectorWithPointers(_list_elements);
    };

    virtual size_t getDimension() const
    {
      return N_DIM;
    };
    virtual const CoordinateSystem& getCoordinateSystem() const
    {
        return _coord_system;        
    };

    virtual double getMinEdgeLength() const 
    {
        std::cout << "***Warning: getMinEdgeLength() just returns 0." << std::endl;
        return .0;
    };

    virtual size_t getNumberOfNodes() const { return _list_nodes.size(); };

    virtual Node* getNode( size_t id ) const {
        return _list_nodes.at(id);
    };

    size_t setNode( size_t node_id, GeoLib::Point &x ) {
        size_t new_node_id = node_id;
        if (node_id<_list_nodes.size()) {
            INode *node = this->_list_nodes.at(node_id);
            node->setX(x);
        } else {
            new_node_id = this->_list_nodes.size();
            this->_list_nodes.push_back(new Node(new_node_id, x));
        }
        return new_node_id;
    };

    virtual void getNodeCoordinates(size_t node_id, GeoLib::Point &pt) const {
        INode* nod = _list_nodes[node_id];
        pt = *nod->getData();
    }

    const std::vector<Node*>& getNodeVector() const { return _list_nodes;};

    virtual size_t getNumberOfElements() const { return _list_elements.size(); };

    virtual IElement* getElemenet( size_t element_id ) const {
        assert(element_id<_list_elements.size());
        return _list_elements.at(element_id);
    };

    size_t addElement( IElement *e) {
        e->setElementID(_list_elements.size());
        _list_elements.push_back(e);
        return e->getElementID();
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

typedef TemplateUnstructuredMesh<1, CoordinateSystem::X> UnstructuredMesh1d;
typedef TemplateUnstructuredMesh<2, CoordinateSystem::XY> UnstructuredMesh2d;
typedef TemplateUnstructuredMesh<3, CoordinateSystem::XYZ> UnstructuredMesh3d;

}
