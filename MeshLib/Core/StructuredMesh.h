
#pragma once

#include "GeoLib/Core/Point.h"
#include "IMesh.h"

namespace MeshLib
{
/**
 * \brief a structured mesh class
 */
class StructuredMesh : public IMesh
{
private:
    GeoLib::Point _origin[3];
    GeoLib::Point _length[3];
    GeoLib::Point  _unit_length[3];
    size_t  _number_of_nodes_per_dimension[3];
    IElement* _e;
    INode* _nod;

public:
    StructuredMesh(CoordinateSystem::CoordinateSystemType coord, GeoLib::Point &org_pt, GeoLib::Point &len, GeoLib::Point &unit_len);
    virtual ~StructuredMesh();

    /// get the number of nodes
    size_t getNumberOfNodes() const {
      return _number_of_nodes_per_dimension[0] + _number_of_nodes_per_dimension[1] + _number_of_nodes_per_dimension[2];
    }

    /// get a node object
    INode* getNode( size_t id ) const {
      return _nod;
    }

    /// get the number of elements
    size_t getNumberOfElements() const {
      return _unit_length[0] + _unit_length[1] + _unit_length[2]; 
    }

    /// get an element
    IElement* getElemenet( size_t element_id ) const {
      //set e
      return _e;
    }
};

}
