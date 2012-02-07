
#pragma once

#include <algorithm>

#include "Base/MemoryTools.h"
#include "GeoLib/Core/Point.h"

#include "IMesh.h"
#include "IElement.h"
#include "ElementFactory.h"
#include "MeshGeometricProperties.h"

namespace MeshLib
{

/**
 * \brief a structured mesh class
 */
template <ElementType::type T_ELEMENT_TYPE>
class StructuredMesh : public IMesh
{
private:
    GeoLib::Point _origin;
    double _length[3];
    double  _unit_length[3];
    size_t  _number_of_elements_per_dimension[3];
    size_t  _number_of_nodes_per_dimension[3];
    IElement* _e;
    INode* _nod;
   size_t _ele_size;
    size_t _nod_size;
    std::vector<IElement*> _list_edge_elements;
    MeshGeometricProperty _geo_prop;
    bool _isAxisymmetric;

    void construct();
    void getNodeCoordinates(size_t id, GeoLib::Point* p) const;
public:

    ///
    StructuredMesh(CoordinateSystemType::type coord, GeoLib::Point &org_pt, const double* len, const double* unit_len) : _origin(org_pt)
    {
        _isAxisymmetric = false;
        _geo_prop.setCoordinateSystem(coord);
        const size_t dim = getDimension();
        for (size_t i=0; i<dim; i++)
            _length[i] = len[i];
        for (size_t i=0; i<dim; i++)
            _unit_length[i] = unit_len[i];

        double d = _unit_length[0];
        for (size_t i=1; i<getDimension(); i++)
            d = std::min(d, _unit_length[i]);
        _geo_prop.setMinEdgeLength(d);

        construct();
    }

    virtual ~StructuredMesh()
    {
        delete _e;
        delete _nod;
        Base::destroyStdVectorWithPointers(_list_edge_elements);
    }

    size_t getDimension() const 
    {
        return _geo_prop.getCoordinateSystem()->getDimension();
    }

    /// get the number of elements
    size_t getNumberOfElements() const {
        return _ele_size;
    }
    /// get an element
    IElement* getElemenet( size_t element_id ) const;

    /// get the number of nodes
    size_t getNumberOfNodes() const {
        return _nod_size;
    }
    /// get a node object
    INode* getNode( size_t id ) const;
    ///
    /// get a point object
    const GeoLib::Point* getNodeCoordinatesRef(size_t id) const 
    {
        getNodeCoordinates(id, const_cast<GeoLib::Point*>(_nod->getData()));
        return _nod->getData();
    }
    /// get a point object
    GeoLib::Point getNodeCoordinates(size_t id) const 
    {
        return *getNodeCoordinatesRef(id);
    }
    void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const 
    {
        for (size_t i=0; i<vec_node_id.size(); i++)
            vec_pt.push_back(*this->getNodeCoordinatesRef(vec_node_id[i]));
    }
    
    /// add a new element
    void addEdgeElement(IElement* e) {
        _list_edge_elements.push_back(e);
    }

    const MeshGeometricProperty* getGeometricProperty() const
    {
        return &_geo_prop;
    }

    void setAxisymmetric(bool flag)
    {
        _isAxisymmetric = flag;
    }
    bool isAxisymmetric() const
    {
        return _isAxisymmetric;
    }

};

template<> void StructuredMesh<ElementType::QUAD>::construct();
template<> IElement* StructuredMesh<ElementType::QUAD>::getElemenet( size_t element_id ) const; 
template<> INode* StructuredMesh<ElementType::QUAD>::getNode( size_t id ) const; 
template<> void StructuredMesh<ElementType::QUAD>::getNodeCoordinates(size_t id,  GeoLib::Point* p) const;

}
