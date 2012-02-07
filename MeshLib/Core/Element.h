
#pragma  once

#include <vector>
#include <algorithm>

#include "Base/MemoryTools.h"
#include "GeoLib/MathTools.h"
#include "IElement.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "IMesh.h"

namespace MeshLib
{

/**
 * 
 */
template <ElementType::type TYPE, size_t NUMBER_OF_NODES, size_t DIMENSION, size_t NUMBER_OF_FACES, size_t NUMER_OF_EDGES>
class TemplateElement : public IElement
{
private:
    size_t _list_node_id[NUMBER_OF_NODES];
    std::vector<IElement*> _list_edges; // memory is assigned later when required

public:
    TemplateElement() {
        initialize();
    };

    TemplateElement(size_t element_id) {
        initialize();
        this->setID(element_id);
    };

    virtual ~TemplateElement() 
    {
    };

    void initialize() {
        for (int i=0; i<NUMBER_OF_NODES; i++)
            _list_node_id[i] = 0;
    }

    virtual bool operator==(IElement& e) 
    {
        throw std::exception("Element::== is not implemented yet.");
        return false;
    }

    virtual bool hasNodeIds(std::vector<size_t> &sorted_node_ids) const
    {
        std::vector<size_t> vec(_list_node_id, _list_node_id+NUMBER_OF_NODES);
        std::sort (vec.begin(), vec.end());
        return std::equal (vec.begin(), vec.end(), sorted_node_ids.begin());
    }

    ElementType::type getElementType() const {return TYPE;};
    size_t getNumberOfNodes() const {return NUMBER_OF_NODES;};
    size_t getDimension() const {return DIMENSION;};
    size_t getNumberOfFaces() const {return NUMBER_OF_FACES;};
    size_t getNumberOfEdges() const {return NUMER_OF_EDGES;};

    IElement* getEdgeElement(size_t edge_id) 
    {
        if (_list_edges.size()<NUMER_OF_EDGES)
            return 0;

        return _list_edges[edge_id];
    }

    void setEdgeElement(size_t edge_id, IElement* e) 
    {
        if (_list_edges.size()==0)
            _list_edges.resize(NUMER_OF_EDGES, 0);
        _list_edges[edge_id] = e;
    }

    virtual void setNodeID(size_t local_node_id, size_t node_id) {
        assert (local_node_id < NUMBER_OF_NODES);
        _list_node_id[local_node_id] = node_id;
    }

    virtual size_t getNodeID(size_t local_node_id) const {
        assert (local_node_id < NUMBER_OF_NODES);
        return _list_node_id[local_node_id];
    };
};

/**
 * \brief Unstructured element
 */
template <ElementType::type TYPE, size_t NUMBER_OF_NODES, size_t DIMENSION, size_t NUMBER_OF_FACES, size_t NUMER_OF_EDGES>
class TemplateUnstructuredElement : public TemplateElement<TYPE, NUMBER_OF_NODES, DIMENSION, NUMBER_OF_FACES, NUMER_OF_EDGES>
{
private:
    //INode* _list_nodes[NUMBER_OF_NODES];
    IElementCoordinatesMapping *_coord_map;

public:
    TemplateUnstructuredElement() {
        initialize();
    };
    TemplateUnstructuredElement(size_t element_id) {
        initialize();
        this->setID(element_id);
    }
    virtual ~TemplateUnstructuredElement() 
    {
    };

    void initialize() {
        _coord_map = 0;
    }

    virtual void getNodeIDsOfEdgeElement(size_t edge_id, std::vector<size_t> &vec_node_ids) const
    {
        vec_node_ids.resize(0);
    }

    virtual ElementType::type getEdgeElementType(size_t edge_id) const 
    {
        return ElementType::INVALID;
    }


    //void setNode(size_t local_node_id, INode* nod) {
    //    _list_nodes[local_node_id] = nod;
    //}

    //INode* getNode(size_t local_node_id) const {
    //    return _list_nodes[local_node_id];
    //}

    //virtual const GeoLib::Point* getNodeCoordinates( size_t i_nod ) const 
    //{
    //    return _list_nodes[i_nod]->getData();
    //};

    virtual void setMappedCoordinates(IElementCoordinatesMapping* mapping) 
    {
        _coord_map = mapping;
    }

    virtual IElementCoordinatesMapping* getMappedCoordinates() 
    {
        return _coord_map;
    };
};

// elements
typedef TemplateUnstructuredElement<ElementType::LINE, 2, 1, 2, 0> Line;
typedef TemplateUnstructuredElement<ElementType::TRIANGLE, 3, 2, 3, 3> TemplateTriangle;
class Triangle : public TemplateTriangle
{
public:
    //double getArea() const
    //{
    //    return GeoLib::triangleArea(*getNodeCoordinates(0), *getNodeCoordinates(1), *getNodeCoordinates(2));
    //}
};
typedef TemplateUnstructuredElement<ElementType::QUAD, 4, 2, 4, 4> Quadrirateral;
typedef TemplateUnstructuredElement<ElementType::TETRAHEDRON, 4, 3, 4, 6> Tetrahedron;
typedef TemplateUnstructuredElement<ElementType::PYRAMID, 5, 3, 5, 8> Pyramid;
typedef TemplateUnstructuredElement<ElementType::PRISM, 6, 3, 5, 9> Prism;
typedef TemplateUnstructuredElement<ElementType::HEXAHEDRON, 8, 3, 6, 12> Hexahedron;

template <> 
void TemplateTriangle::getNodeIDsOfEdgeElement(size_t edge_id, std::vector<size_t> &vec_node_ids) const 
{
    vec_node_ids.resize(2);
    vec_node_ids[0] = this->getNodeID(edge_id);
    vec_node_ids[1] = this->getNodeID((edge_id+1)%3);
}

template <> 
void Quadrirateral::getNodeIDsOfEdgeElement(size_t edge_id, std::vector<size_t> &vec_node_ids) const 
{
    vec_node_ids.resize(2);
    vec_node_ids[0] = this->getNodeID(edge_id);
    vec_node_ids[1] = this->getNodeID((edge_id+1)%4);
}

template <> 
ElementType::type TemplateTriangle::getEdgeElementType(size_t edge_id) const
{
    return ElementType::LINE;
}

template <> 
ElementType::type Quadrirateral::getEdgeElementType(size_t edge_id) const
{
    return ElementType::LINE;
}

///**
// * \brief Analytical element
// */
//class AnalyticalQuadrirateral : public Quadrirateral
//{
//private:
//    IMesh *_msh;
//public:
//    AnalyticalQuadrirateral(IMesh *msh) 
//    {
//        _msh = msh;
//    }
//
//    const GeoLib::Point* getNodeCoordinates( size_t i_nod ) const 
//    {
//        return _msh->getNode(i_nod)->getData();
//    };
//};


} // end namespace

