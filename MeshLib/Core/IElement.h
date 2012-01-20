
#pragma  once

#include <vector>

namespace MeshLib
{

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
    void setGroupID(size_t id) {_group_id = id;};
    size_t getGroupID() const {return _group_id;};

    virtual ElementType::type getElementType() const = 0;
    virtual size_t getElementDimension() const = 0;
    virtual size_t getNumberOfNodes() const = 0;
    virtual size_t getNumberOfFaces() const = 0;
    virtual size_t getNumberOfEdges() const = 0;

    virtual void setNodeID(size_t local_node_id, size_t node_id) = 0;
    virtual size_t getNodeID(size_t local_node_id) const = 0;

    void getNodeIDList( std::vector<size_t> &e_node_id_list ) const
    {
        e_node_id_list.resize(this->getNumberOfNodes());
        for (int i=0; i<this->getNumberOfNodes(); i++)
            e_node_id_list[i] = this->getNodeID(i);
    };

    double getLocalX( size_t i_nod ) 
    {
        throw std::exception("The method or operation is not implemented.");
    }



};
}
