
#pragma once

#include <vector>

#include "GeoLib/Core/Point.h"
#include "IElement.h"

namespace MeshLib
{
class CoordinateSystem;
class INode;
class MeshGeometricProperty;

/**
 * \brief Interface to mesh classes
 *
 * Mesh classes should have the following data,
 * - Node
 * - Vertex
 * - Edge
 * - Face
 * - Element
 */
class IMesh
{
public:
    virtual void setID(size_t id) = 0;
    virtual size_t getID() const = 0;
    ///// get coordinate systems
    //virtual const CoordinateSystem* getCoordinateSystem() const = 0;
    ///// set coordinate systems
    //virtual void setCoordinateSystem(CoordinateSystem &coord) = 0;
    virtual size_t getDimension() const = 0;
    /// get mesh geometric property
    virtual const MeshGeometricProperty* getGeometricProperty() const = 0;
    /// return if this mesh is axisymmetric or not
    virtual bool isAxisymmetric() const = 0;
    virtual void setAxisymmetric(bool flag) = 0;

    /// get the number of elements
    virtual size_t getNumberOfElements() const = 0;
    /// get the number of elements of the given type
    virtual size_t getNumberOfElements(ElementShape::type ele_type) const = 0;
    /// get an element
    virtual IElement* getElemenet( size_t element_id ) const = 0;

    /// get the number of nodes
    virtual size_t getNumberOfNodes() const = 0;
    /// get a node object
    virtual INode* getNode( size_t id ) const = 0;
    /// get a point object
    virtual const GeoLib::Point* getNodeCoordinatesRef(size_t id) const = 0;
    virtual GeoLib::Point getNodeCoordinates(size_t id) const = 0;
    /// get a list of points in the given element
    virtual void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const = 0;


    /// add a new element
    virtual void addEdgeElement(IElement*) = 0;
    virtual size_t getNumberOfEdges() const = 0;
    virtual IElement* getEdgeElement(size_t edge_id) = 0;
};


class IMixedOrderMesh : public IMesh
{
public:
    //#
    virtual size_t getMaxiumOrder() const = 0;
    virtual void setCurrentOrder(size_t order) = 0;
    virtual size_t getCurrentOrder() const = 0;
    virtual size_t getNumberOfNodes() const = 0; // i don't understand why this redefinition is needed
    virtual size_t getNumberOfNodes(size_t order) const = 0;
    size_t getNumberOfTotalNodes() const
    {
        return getNumberOfNodes(getMaxiumOrder());
    };
};

}
