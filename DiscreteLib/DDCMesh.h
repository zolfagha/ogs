
#pragma once

#include "DomainDecomposition.h"
#include "MeshLib/Core/IMesh.h"

namespace DiscreteLib
{

class AbstractDDCMesh : public MeshLib::IMesh
{
protected:
    DDCGlobal* _ddc;
    std::vector<MeshLib::IMesh*> _list_local_msh;
    MeshLib::IMesh* _msh0;
public:
    AbstractDDCMesh(DDCGlobal& ddc) : _ddc(&ddc)
    {
        for (size_t i=0; i<ddc.getNumberOfSubDomains(); i++)
            _list_local_msh.push_back(ddc.getSubDomain(i)->getLoalMesh());
        _msh0 = _list_local_msh[0];
    }
    virtual ~AbstractDDCMesh() {};

    // ok
    size_t virtual getDimension() const { return _msh0->getDimension(); }
    bool isAxisymmetric() const {return _msh0->isAxisymmetric();}
    void setAxisymmetric(bool flag) {return _msh0->setAxisymmetric(flag);}

    // need to implement
    const MeshLib::MeshGeometricProperty* getGeometricProperty() const
    {
        return _msh0->getGeometricProperty();
    }


    /// not supporting methods
    void addEdgeElement( MeshLib::IElement*) {};
    size_t getNumberOfEdges() const {return 0;};
    MeshLib::IElement* getEdgeElement(size_t edge_id) {return 0;};
};


class NodeDDCMesh : public AbstractDDCMesh
{
public:
    NodeDDCMesh(DDCGlobal& ddc) : AbstractDDCMesh(ddc)
    {
    }
    virtual ~NodeDDCMesh() {};

    size_t getNumberOfElements() const 
    {
        size_t cnt = 0;
        for (size_t i=0; i<_ddc->getNumberOfSubDomains(); i++)
        {
            DDCSubDomain* dom = _ddc->getSubDomain(i);
        }
        return cnt;
    }

    size_t getNumberOfElements( MeshLib::ElementShape::type ele_type) const
    {
        size_t cnt = 0;
        return cnt;
    }

    MeshLib::IElement* getElemenet( size_t element_id ) const
    {

    }

    size_t getNumberOfNodes() const
    {
        size_t cnt = 0;
        for (size_t i=0; i<_list_local_msh.size(); i++)
            cnt += _list_local_msh[i]->getNumberOfNodes();
        return cnt;
    }

    MeshLib::INode* getNode( size_t id ) const
    {
        size_t domid = _ddc->findDomainID(id);
        return _list_local_msh[domid]->getNode(node_address(id, domid));
    }

    const GeoLib::Point* getNodeCoordinatesRef(size_t id) const
    {
        size_t domid = _ddc->findDomainID(id);
        return _list_local_msh[domid]->getNodeCoordinatesRef(node_address(id, domid));
    }
    GeoLib::Point getNodeCoordinates(size_t id) const
    {
        size_t domid = _ddc->findDomainID(id);
        return _list_local_msh[domid]->getNodeCoordinates(node_address(id, domid));
    }
    void getListOfNodeCoordinates(const std::vector<size_t> &vec_node_id, std::vector<GeoLib::Point> &vec_pt) const
    {
        //TODO implement
    }

private:
    size_t node_address(size_t global_id, size_t domid) const
    {
        size_t local_id = _ddc->getSubDomain(domid)->getGlobalLocalIdMap()->global2local(global_id);
        return local_id;
    };
};

class ElementDDCMesh : public AbstractDDCMesh
{
public:
    /// get the number of elements
    size_t getNumberOfElements() const 
    {
        size_t cnt = 0;
        for (size_t i=0; i<_list_local_msh.size(); i++)
            cnt += _list_local_msh[i]->getNumberOfElements();
        return cnt;
    }

    /// get the number of elements of the given type
    size_t getNumberOfElements( MeshLib::ElementShape::type ele_type) const
    {
        size_t cnt = 0;
        for (size_t i=0; i<_list_local_msh.size(); i++)
            cnt += _list_local_msh[i]->getNumberOfElements(ele_type);
        return cnt;
    }
};

} //END
