
#pragma once

#include <vector>

#include "GeoLib/Core/Point.h"


namespace MeshLib
{
class CoordinateSystem;
class INode;
class IElement;
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
    ///// get coordinate systems
    //virtual const CoordinateSystem* getCoordinateSystem() const = 0;
    ///// set coordinate systems
    //virtual void setCoordinateSystem(CoordinateSystem &coord) = 0;
    size_t virtual getDimension() const = 0;

    /// get the number of elements
    virtual size_t getNumberOfElements() const = 0;
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

    /// get mesh geometric property
    virtual const MeshGeometricProperty* getGeometricProperty() const = 0;

    /// return if this mesh is axisymmetric or not
    virtual bool isAxisymmetric() const = 0;
    virtual void setAxisymmetric(bool flag) = 0;
};

}
