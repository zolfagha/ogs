
#pragma once

namespace MeshLib
{

class IElement;
class Node;
    
class IMesh
{
public:
    virtual ~IMesh(){};

    virtual size_t getNumberOfNodes() const = 0;
    virtual size_t getNumberOfElements() const = 0;
    virtual void getNodeCoordinates( size_t node_id, double pt[3] ) const = 0;
    virtual IElement* getElemenet( size_t element_id ) const = 0;

    virtual Node* getNode( size_t id ) const = 0;
};
}
