
#pragma once

#include "INode.h"
#include "IElement.h"
#include "CoordinateSystem.h"

namespace MeshLib
{

/**
 * \brief Interface to mesh classes
 *
 * 
 */
class IMesh
{
public:
    /// get coordinate systems
    virtual const CoordinateSystem* getCoordinateSystem() const = 0;
    /// set coordinate systems
    virtual void setCoordinateSystem(CoordinateSystem &coord) = 0;

    /// get the number of elements
    virtual size_t getNumberOfElements() const = 0;
    /// get an element
    virtual IElement* getElemenet( size_t element_id ) const = 0;

    /// get the number of nodes
    virtual size_t getNumberOfNodes() const = 0;
    /// get a node object
    virtual INode* getNode( size_t id ) const = 0;


    /// add a new element
    virtual size_t addEdgeElement(IElement*) = 0;
};

}
