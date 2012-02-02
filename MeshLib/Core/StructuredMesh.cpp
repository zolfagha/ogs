
#include "StructuredMesh.h"

namespace MeshLib
{

template<> void StructuredMesh<ElementType::QUAD>::construct()
{
    const size_t dim = this->_coord.getDimension();
    for (size_t i=0; i<dim; i++)
        _number_of_elements_per_dimension[i] = static_cast<size_t>(_length[i] / _unit_length[i]);

    _ele_size = 1;
    for (size_t i=0; i<dim; i++)
        _ele_size *= _number_of_elements_per_dimension[i];

    for (size_t i=0; i<dim; i++)
        _number_of_nodes_per_dimension[i] = _number_of_elements_per_dimension[i] + 1;

    _nod_size = 1;
    for (size_t i=0; i<dim; i++)
        _nod_size *= _number_of_nodes_per_dimension[i];

    _e = new Quadrirateral();
    _nod = new Node();
}


template<> IElement* StructuredMesh<ElementType::QUAD>::getElemenet( size_t element_id ) const 
{
    //set e
    _e->initialize();
    _e->setID(element_id);
    const size_t x_j = element_id / _number_of_elements_per_dimension[0];
    const size_t offset_y1 = x_j*_number_of_nodes_per_dimension[0];
    const size_t offset_y2 = (x_j+1)*_number_of_nodes_per_dimension[0];
    const size_t k = element_id % _number_of_elements_per_dimension[0];
    for (size_t i=0; i<4; i++) {
        _e->setNodeID(0, offset_y1+k);
        _e->setNodeID(1, offset_y1+k+1);
        _e->setNodeID(2, offset_y2+k+1);
        _e->setNodeID(3, offset_y2+k);
        //for (size_t l=0; l<4; l++)
        //    e->setNode(l, msh->getNode(e->getNodeID(l)));
    }
    return _e;
};

template<> INode* StructuredMesh<ElementType::QUAD>::getNode( size_t id ) const 
{
    _nod->setNodeID(id);
    getNodeCoordinates(id, const_cast<GeoLib::Point*>(_nod->getData()));

    return _nod;
};

template<> void StructuredMesh<ElementType::QUAD>::getNodeCoordinates(size_t id,  GeoLib::Point* p) const
{
    double *pt = const_cast<double*>(p->getData());
    size_t k_x = id % _number_of_nodes_per_dimension[0]; 
    size_t j_y = id / _number_of_nodes_per_dimension[1]; 
    pt[0] = _unit_length[0]*k_x + _origin[0];
    pt[1] = _unit_length[1]*j_y + _origin[1];
    pt[2] = _origin[2];
}

}
