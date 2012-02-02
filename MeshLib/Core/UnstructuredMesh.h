
#pragma once

#include <iostream>
#include <vector>
#include <limits>

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
    std::vector<IElement*> _list_edge_elements;
    CoordinateSystem _coord_system;

public:
    TemplateUnstructuredMesh() : _coord_system(coord_type)
    {
    };
    virtual ~TemplateUnstructuredMesh(){
        destroyStdVectorWithPointers(_list_nodes);
        destroyStdVectorWithPointers(_list_elements);
        destroyStdVectorWithPointers(_list_edge_elements);
    };

    virtual size_t getDimension() const
    {
      return N_DIM;
    };
    virtual const CoordinateSystem* getCoordinateSystem() const
    {
        return &_coord_system;        
    };
    /// set coordinate systems
    virtual void setCoordinateSystem(CoordinateSystem &coord) 
    {
        _coord_system = coord;
    }

    virtual double getMinEdgeLength() const 
    {
        std::cout << "***Warning: getMinEdgeLength() just returns numeric_limits<double>::epsilon()." << std::endl;
        return std::numeric_limits<double>::epsilon();
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

    const GeoLib::Point*  getNodeCoordinatesRef(size_t id) const
    {
        return getNode(id)->getData();
    }

    GeoLib::Point getNodeCoordinates(size_t id) const
    {
        return *getNode(id)->getData();
    }

    virtual void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const 
    {
        for (size_t i=0; i<vec_node_id.size(); i++)
            vec_pt.push_back(*this->getNodeCoordinatesRef(vec_node_id[i]));
    }

    const std::vector<Node*>& getNodeVector() const { return _list_nodes;};

    virtual size_t getNumberOfElements() const { return _list_elements.size(); };

    virtual IElement* getElemenet( size_t element_id ) const {
        assert(element_id<_list_elements.size());
        return _list_elements.at(element_id);
    };

    size_t addElement( IElement *e) {
        e->setID(_list_elements.size());
        _list_elements.push_back(e);
        return e->getID();
    };

    void addEdgeElement(IElement *e) 
    {
        e->setID(_list_edge_elements.size());
        _list_edge_elements.push_back(e);
    }


    void construct() {

        //// set 
        //for (size_t i=0; i<getNumberOfElements(); i++) {
        //    IElement *e = getElemenet(i);
        //    for (size_t j=0; j<e->getNumberOfNodes; j++)
        //        e->setNode(j, getNode(e->getNodeID(j)));
        //}

        ////set node connectivity
        //const size_t n_ele = this->getNumberOfElements();
        //for (size_t i=0; i<n_ele; i++) {
        //    const IElement* e = _list_elements[i];
        //    for (size_t j=0; j<e->getNumberOfNodes(); j++) {
        //        Node* nod_j = getNode(e->getNodeID(j));
        //        for (size_t k=j+1; k<e->getNumberOfNodes(); k++) {
        //            Node* nod_k = getNode(e->getNodeID(k));
        //            nod_j->addConnectedNode(e->getNodeID(k));
        //            nod_k->addConnectedNode(e->getNodeID(j));
        //        }
        //    };
        //}
    };
};

typedef TemplateUnstructuredMesh<1, CoordinateSystem::X> UnstructuredMesh1d;
typedef TemplateUnstructuredMesh<2, CoordinateSystem::XY> UnstructuredMesh2d;
typedef TemplateUnstructuredMesh<3, CoordinateSystem::XYZ> UnstructuredMesh3d;


template<size_t N_DIM, CoordinateSystem::CoordinateSystemType coord_type>
class TemplateMixedOrderUnstructuredMesh : public TemplateUnstructuredMesh<N_DIM, coord_type>
{
private:
    size_t _order;
public:
    TemplateMixedOrderUnstructuredMesh() : _order(1) {};

    void setOrder(size_t order)
    {
        _order = order;
    };

    size_t getOrder() const
    {
        return _order;
    };
};

}
