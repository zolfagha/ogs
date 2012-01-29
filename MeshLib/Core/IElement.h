
#pragma  once

#include <vector>

#include "GeoLib/Core/Point.h"

#include "MeshLib/Core/INode.h"
#include "CoordinateSystem.h"

namespace MeshLib
{

/// List of available element types
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

class IElementCoordinatesMapping;

/**
 * Interface of element classes
 */
class IElement
{
private:
    size_t _element_id;
    size_t _group_id;
public:
    IElement():_element_id(0), _group_id(0) {};
    virtual ~IElement() {};

    virtual bool operator==(IElement& e) = 0;
    virtual bool hasNodeIds(std::vector<size_t> &node_ids) const = 0;

    size_t getID() const {return _element_id;};
    void setID(size_t id) {_element_id = id;};

    size_t getGroupID() const {return _group_id;};
    void setGroupID(size_t id) {_group_id = id;};

    virtual ElementType::type getElementType() const = 0;

    virtual size_t getDimension() const = 0;

    virtual size_t getNumberOfNodes() const = 0;
    virtual void setNodeID(size_t local_node_id, size_t node_id) = 0;
    virtual size_t getNodeID(size_t local_node_id) const = 0;
    virtual void setNode(size_t local_node_id, INode*) = 0;
    virtual INode* getNode(size_t local_node_id) const = 0;
    void getNodeIDList( std::vector<size_t> &e_node_id_list ) const
    {
        e_node_id_list.resize(this->getNumberOfNodes());
        for (size_t i=0; i<this->getNumberOfNodes(); i++)
            e_node_id_list[i] = this->getNodeID(i);
    };
    virtual const GeoLib::Point* getNodeCoordinates( size_t i_nod ) const = 0;

    virtual size_t getNumberOfFaces() const = 0;

    virtual size_t getNumberOfEdges() const = 0;
    virtual IElement* getEdgeElement(size_t edge_id) = 0;
    virtual void setEdgeElement(size_t edge_id, IElement* e) = 0;
    virtual void getNodeIDsOfEdgeElement(size_t edge_id, std::vector<size_t> &vec_node_ids) const = 0;
    virtual ElementType::type getEdgeElementType(size_t edge_id) const = 0;

    virtual void setMappedCoordinates(IElementCoordinatesMapping* mapping) = 0;
    virtual IElementCoordinatesMapping* getMappedCoordinates() = 0;

};

}
