
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
#include "MeshGeometricProperties.h"


namespace MeshLib
{

/**
 * \brief Unstructured mesh class
 *
 *
 */
class UnstructuredMesh : public IMesh
{
private:
    std::vector<Node*> _list_nodes;
    std::vector<IElement*> _list_elements;
    std::vector<IElement*> _list_edge_elements;
    MeshGeometricProperty _geo_prop;
    bool _isAxisymmetric;

public:
    UnstructuredMesh() 
    {
    };

    UnstructuredMesh(const CoordinateSystemType::type coord) 
    {
        _geo_prop.setCoordinateSystem(coord);
    };

    virtual ~UnstructuredMesh()
    {
        Base::destroyStdVectorWithPointers(_list_nodes);
        Base::destroyStdVectorWithPointers(_list_elements);
        Base::destroyStdVectorWithPointers(_list_edge_elements);
    }

    //------------------------------------------------------------------------
    size_t getDimension() const
    {
        return _geo_prop.getCoordinateSystem()->getDimension();
    };

    /// get mesh geometric property
    const MeshGeometricProperty* getGeometricProperty() const
    {
        return &_geo_prop;
    }

    MeshGeometricProperty* getGeometricProperty()
    {
        return &_geo_prop;
    }

    bool isAxisymmetric() const
    {
        return _isAxisymmetric;
    }
    void setAxisymmetric(bool flag)
    {
        _isAxisymmetric = flag;
    }


    //------------------------------------------------------------------------
    //
    virtual size_t getNumberOfNodes() const 
    { 
        return _list_nodes.size(); 
    };

    virtual Node* getNode( size_t id ) const 
    {
        return _list_nodes.at(id);
    };

    size_t addNode( GeoLib::Point &x ) 
    {
        size_t new_node_id = _list_nodes.size();
        this->_list_nodes.push_back(new Node(new_node_id, x));
        return  new_node_id;
    }

    size_t setNode( size_t node_id, GeoLib::Point &x ) 
    {
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

    const std::vector<Node*>& getNodeVector() const 
    { 
        return _list_nodes;
    };

    //------------------------------------------------------------------------
    virtual size_t getNumberOfElements() const 
    { 
        return _list_elements.size(); 
    };

    virtual IElement* getElemenet( size_t element_id ) const 
    {
        assert(element_id<_list_elements.size());
        return _list_elements.at(element_id);
    };

    size_t addElement( IElement *e) 
    {
        e->setID(_list_elements.size());
        _list_elements.push_back(e);
        return e->getID();
    };

    void addEdgeElement(IElement *e) 
    {
        e->setID(_list_edge_elements.size());
        _list_edge_elements.push_back(e);
    }
};


/**
 * \brief Mixed order mesh
 */
class MixedOrderUnstructuredMesh : public UnstructuredMesh
{
private:
    size_t _order;
    size_t _max_order;
    std::vector<size_t> _n_nodes;


public:
    MixedOrderUnstructuredMesh() : _order(1) {};


    void setCurrentOrder(size_t order)
    {
        _order = order;
    };

    size_t getCurrentOrder() const
    {
        return _order;
    };

    size_t getMaxiumOrder() const
    {
        return _max_order;
    }

    /// get an element
    virtual IElement* getElemenet( size_t element_id ) const = 0;

    /// get the number of nodes
    size_t getNumberOfNodes() const
    {
        return getNumberOfNodes(_order);
    }
    size_t getNumberOfNodes(size_t order) const
    {
        return _n_nodes[order];
    }

};

class MixedOrderElement
{
public:
    size_t getNumberOfNodes(size_t order);
};

}
