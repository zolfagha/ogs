
#pragma once

#include "INode.h"
#include "IElement.h"

namespace MeshLib
{

template<typename Tpos>
class IMesh
{
public:
    virtual ~IMesh(){};

    virtual size_t getNumberOfNodes() const = 0;
    virtual size_t getNumberOfElements() const = 0;
    virtual void getNodeCoordinates( size_t node_id, Tpos &x ) const = 0;
    virtual IElement* getElemenet( size_t element_id ) const = 0;

    virtual INode<Tpos>* getNode( size_t id ) const = 0;
};
}
