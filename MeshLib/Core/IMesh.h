
#pragma once

#include "MathLib/Vector.h"

#include "INode.h"
#include "IElement.h"

namespace MeshLib
{

class IMesh
{
public:
    virtual ~IMesh(){};

    virtual size_t getDimension() const = 0;

    virtual size_t getNumberOfNodes() const = 0;
    virtual INode* getNode( size_t id ) const = 0;

    virtual size_t getNumberOfElements() const = 0;
    virtual IElement* getElemenet( size_t element_id ) const = 0;

};

}
