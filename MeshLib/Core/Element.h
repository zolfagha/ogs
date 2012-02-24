
#pragma  once

#include <vector>
#include <algorithm>

#include "Base/CodingTools.h"
#include "GeoLib/MathTools.h"

#include "IElement.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "IMesh.h"
#include "ElementTopology.h"

namespace MeshLib
{

/**
 * 
 */
template <typename T_TOPO>
class TemplateElement : public IElement
{
public:
    TemplateElement() 
    {
        initialize();
    };

    TemplateElement(size_t element_id) 
    {
        initialize();
        this->setID(element_id);
    };

    virtual ~TemplateElement() 
    {
    };

    void initialize() 
    {
        _element_id = 0;
        _group_id = 0;
        _order = 1;
        _max_order = _order;
        _coord_map = 0;
        _list_node_id.resize(getNumberOfNodes(_order), 0);
        std::fill(_list_node_id.begin(), _list_node_id.end(), 0);
        _list_edges.clear();
        _list_edges.resize(getNumberOfEdges(), 0);
    }

    IElement* clone() const
    {
        TemplateElement<T_TOPO>* obj = new TemplateElement<T_TOPO>();
        obj->_group_id = _group_id;
        obj->_order = _order;
        obj->_list_node_id.assign(_list_node_id.begin(), _list_node_id.end());
        //obj->_coord_map = _coord_map;
        return obj;
    }

    /// return this element id
    size_t getID() const {return _element_id;};
    /// set this element id
    void setID(size_t id) {_element_id = id;};

    /// return the group id of this element
    size_t getGroupID() const {return _group_id;};
    /// set the group if of this element
    void setGroupID(size_t id) {_group_id = id;};

    void setMaximumOrder(size_t order)
    {
        _max_order = order;
        _list_node_id.resize(getNumberOfNodes(order));
    }

    size_t getMaximumOrder() const {return _max_order;};

    /// set current order
    void setCurrentOrder(size_t order) { _order = order; };
    /// get current order
    size_t getCurrentOrder() const { return _order; };

    //virtual bool operator==(IElement& e) 
    //{
    //    if (this->getShapeType()!=e.getShapeType()) return false;
    //    if (this->getNumberOfNodes()!=e.getNumberOfNodes()) return false;
    //    return true;
    //}

    bool hasNodeIds(const std::vector<size_t> &sorted_node_ids) const
    {
        std::vector<size_t> vec(_list_node_id.begin(), _list_node_id.end());
        std::sort (vec.begin(), vec.end());
        return std::equal (vec.begin(), vec.end(), sorted_node_ids.begin());
    }

    ElementShape::type getShapeType() const {return T_TOPO::getShapeType();};
    size_t getDimension() const {return T_TOPO::getDimension();};
    size_t getNumberOfFaces() const {return T_TOPO::getNumberOfFaces();};
    size_t getNumberOfEdges() const {return T_TOPO::getNumberOfEdges();};


    size_t getNumberOfNodes() const 
    {
        return getNumberOfNodes(getCurrentOrder());
    }
    size_t getNumberOfNodes(size_t order) const
    {
        return T_TOPO::getNumberOfNodes(order);
    }

    ElementShape::type getEdgeElementType(size_t edge_id) const
    {
        return T_TOPO::getEdgeElementType(edge_id);
    }

    IElement* getEdgeElement(size_t edge_id) const
    {
        if (_list_edges.size()<getNumberOfEdges())
            return 0;

        return _list_edges[edge_id];
    }

    void setEdgeElement(size_t edge_id, IElement* e) 
    {
        if (_list_edges.size()==0)
            _list_edges.resize(getNumberOfEdges(), 0);
        _list_edges[edge_id] = e;
    }

    void setNodeID(const size_t &local_node_id, const size_t &node_id) 
    {
        //assert (local_node_id < NUMBER_OF_NODES);
        _list_node_id[local_node_id] = node_id;
    }

    size_t getNodeID(size_t local_node_id) const 
    {
        //assert (local_node_id < NUMBER_OF_NODES);
        return _list_node_id[local_node_id];
    };

    /// get a list of node ids
    void getNodeIDList( std::vector<size_t> &e_node_id_list ) const
    {
        getNodeIDList(getCurrentOrder(), e_node_id_list);
    };
    void getNodeIDList( size_t order, std::vector<size_t> &e_node_id_list ) const
    {
        const size_t nnodes = this->getNumberOfNodes(order);
        e_node_id_list.resize(nnodes);
        for (size_t i=0; i<nnodes; i++)
            e_node_id_list[i] = this->getNodeID(i);
    }
    void getListOfNumberOfNodesForAllOrders(std::vector<size_t> &vec) const
    {
        vec.resize(getMaximumOrder());
        for (size_t i=0; i<getMaximumOrder(); i++) {
            vec[i] = this->getNumberOfNodes(i+1);
        }
    }

    void getNodeIDsOfEdgeElement(size_t edge_id, std::vector<size_t> &vec_node_ids) const 
    {
        getNodeIDsOfEdgeElement(getCurrentOrder(), edge_id, vec_node_ids);
    }
    void getNodeIDsOfEdgeElement(size_t order, size_t edge_id, std::vector<size_t> &vec_node_ids) const 
    {
        std::vector<size_t> vec_local_node_ids;
        T_TOPO::getLocalNodeIDsOfEdgeElement(order, edge_id, vec_local_node_ids);
        vec_node_ids.resize(vec_local_node_ids.size());
        for (size_t i=0; i<vec_local_node_ids.size(); i++)
            vec_node_ids[i] = this->getNodeID(vec_local_node_ids[i]);
    };

    void setMappedCoordinates(IElementCoordinatesMapping* mapping) 
    {
        _coord_map = mapping;
    }

    IElementCoordinatesMapping* getMappedCoordinates() 
    {
        return _coord_map;
    };

private:
    size_t _element_id;
    size_t _order;
    size_t _max_order;
    size_t _group_id;
    std::vector<size_t> _list_node_id;
    std::vector<IElement*> _list_edges;
    IElementCoordinatesMapping *_coord_map;
};

// elements
typedef TemplateElement<LineTopology> Line;
typedef TemplateElement<TriangleTopology> Triangle;
typedef TemplateElement<Quad9Topology> Quadrirateral;
typedef TemplateElement<TetraTopology> Tetrahedron;
typedef TemplateElement<PyramidTopology> Pyramid;
typedef TemplateElement<PrismTopology> Prism;
typedef TemplateElement<HexTopology> Hexahedron;


} // end namespace

