
#pragma once

#include <vector>
#include <map>

#include "IMesh.h"
#include "IElement.h"
#include "Node.h"
#include "MeshGeometricProperties.h"


namespace MeshLib
{

class INode;

/**
 * \brief Unstructured mesh class
 *
 *
 */
class UnstructuredMesh : public IMixedOrderMesh
{
private:
    std::vector<Node*> _list_nodes;
    std::vector<IElement*> _list_elements;
    std::map<ElementShape::type, size_t> _n_eles_type;
    std::vector<IElement*> _list_edge_elements;
    MeshGeometricProperty _geo_prop;
    bool _isAxisymmetric;
    size_t _order;
    std::map<size_t, size_t> _map_order_nnodes;

public:
    UnstructuredMesh() : _order(1)
    {
    };

    UnstructuredMesh(const CoordinateSystemType::type coord)  : _order(1)
    {
        _geo_prop.setCoordinateSystem(coord);
    };

    virtual ~UnstructuredMesh();

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
        return _map_order_nnodes.size();
    }

    size_t getNumberOfNodes(size_t order) const
    {
        if (order-1 < _map_order_nnodes.size()) {
            size_t sum = 0;
            for (size_t i=0; i<order; ++i) {
                std::map<size_t, size_t>::const_iterator itr = _map_order_nnodes.find(i+1);
                sum += itr->second;
            }
            return sum;
        }
        return 0;
    }

    //------------------------------------------------------------------------
    size_t getNumberOfNodes() const
    {
        return getNumberOfNodes(_order);
    }

    const std::vector<Node*>& getNodeVector() const 
    { 
        return _list_nodes;
    };

    Node* getNode( size_t id ) const 
    {
        return _list_nodes[id];
    };

    size_t addNode( GeoLib::Point &x, size_t order=1 ); 

    size_t setNodeCoordinates( size_t node_id, GeoLib::Point &x );

    const GeoLib::Point*  getNodeCoordinatesRef(size_t id) const
    {
        return getNode(id)->getData();
    }

    GeoLib::Point getNodeCoordinates(size_t id) const
    {
        return *getNode(id)->getData();
    }

    void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const;


    //------------------------------------------------------------------------
    size_t getNumberOfElements() const 
    { 
        return _list_elements.size(); 
    };

    size_t getNumberOfElements(ElementShape::type ele_type) const;

    IElement* getElemenet( size_t element_id ) const 
    {
        assert(element_id<_list_elements.size());
        IElement* e = _list_elements[element_id];
        e->setCurrentOrder(_order);
        return e;
    };

    size_t addElement( IElement *e);

    //------------------------------------------------------------------------
    void addEdgeElement(IElement *e);
    size_t getNumberOfEdges() const;
    IElement* getEdgeElement(size_t edge_id);
};

}
