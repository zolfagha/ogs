
#pragma once

#include <vector>

#include "Base/CodingTools.h"
#include "Base/BidirectionalMap.h"

#include "MeshLib/Core/IMesh.h"


namespace DiscreteLib
{

struct DecompositionType 
{
    enum type 
    {
        Node,
        Element
    };
};

class DDCSlaveDomain
{
public:
    size_t getDomainID() const {return _dom_id;};
    size_t getNumberOfSlaveObjects() const {return _list_slave_objects.size();};
    size_t getSlaveObjectID(size_t i) const {return _list_slave_objects[i];};
private:
    size_t _dom_id;
    std::vector<size_t> _list_slave_objects;
};

class IDDCGlobaLocalMapping
{
public:
	virtual ~IDDCGlobaLocalMapping() {};
    virtual bool hasGlobal(size_t global_id) = 0;
    virtual size_t global2local(size_t global) = 0;
    virtual size_t local2global(size_t local) = 0;
};

class DDCGlobaLocalMappingOffset : public IDDCGlobaLocalMapping
{
public:
    DDCGlobaLocalMappingOffset(size_t i_start, size_t i_end, size_t offset) : _offset(offset)
    {
        _i_start = i_start;
        _i_end = i_end;
    }
    virtual ~DDCGlobaLocalMappingOffset()
    {
    }
    bool hasGlobal(size_t global_id) 
    {
        return _i_start<=global_id && global_id<_i_end;
    }
    size_t global2local(size_t global)
    {
        return global - _offset;
    }
    size_t local2global(size_t local)
    {
        return local + _offset;
    }
private:
    size_t _i_start;
    size_t _i_end;
    size_t _offset;
};

class DDCGlobaLocalMappingAll : public IDDCGlobaLocalMapping
{
public:
    DDCGlobaLocalMappingAll(Base::BidirectionalMap<size_t, size_t> &map) : _map_global_local(&map)
    {
    }
    virtual ~DDCGlobaLocalMappingAll()
    {
        Base::releaseObject(_map_global_local);
    }
    bool hasGlobal(size_t global_id) 
    {
        return _map_global_local->countInA(global_id)>0;
    }
    size_t global2local(size_t global)
    {
        return _map_global_local->mapAtoB(global);
    }
    size_t local2global(size_t local)
    {
        return _map_global_local->mapBtoA(local);
    }
private:
    Base::BidirectionalMap<size_t, size_t>* _map_global_local;
};

class DDCSubDomain
{
public:
    DDCSubDomain(MeshLib::IMesh &msh, IDDCGlobaLocalMapping &mapping, std::set<size_t>* list_ghosts=0) : _local_msh(&msh), _map_global_local(&mapping)
    {
        _dom_id = 0;
        if (list_ghosts) {
            _list_ghosts.insert(list_ghosts->begin(), list_ghosts->end());
        }
    }
    virtual ~DDCSubDomain()
    {
        Base::releaseObjectsInStdVector(_list_slaves);
    }

    void setDomainID(size_t i) {_dom_id = i;};
    size_t getDomainID() const {return _dom_id;};

    MeshLib::IMesh* getLoalMesh() {return _local_msh;};

    IDDCGlobaLocalMapping* getGlobalLocalIdMap() {return _map_global_local;};

    size_t getNumberOfGhosts() const {return _list_ghosts.size();};
    std::set<size_t>* getGhostList() {return &_list_ghosts;};

    size_t getNumberOfSlaves() const {return _list_slaves.size();};
    DDCSlaveDomain* getSlaveDomain(size_t i) {return _list_slaves[i];};
private:
    size_t _dom_id;
    MeshLib::IMesh* _local_msh;
    IDDCGlobaLocalMapping* _map_global_local;
    std::vector<DDCSlaveDomain*> _list_slaves;
    std::set<size_t> _list_ghosts;
};

class DDCGlobal
{
public:
    DDCGlobal(DecompositionType::type ddc_type) : _ddc_type(ddc_type), _id(0)
    {
        _n_discrete_pt = 0;
    };
    virtual ~DDCGlobal()
    {
        Base::releaseObjectsInStdVector(_list_dom);
    }

    size_t getID() const {return _id;};
    void setID(size_t i) {_id = i;};

    DecompositionType::type getDecompositionType() const {return _ddc_type;};

    size_t addSubDomain(DDCSubDomain* sub)
    {
        _list_dom.push_back(sub);
        sub->setDomainID(_list_dom.size()-1);
        if (_ddc_type==DecompositionType::Node) {
            _n_discrete_pt += sub->getLoalMesh()->getNumberOfNodes();
        } else {
            _n_discrete_pt += sub->getLoalMesh()->getNumberOfElements();
        }
        _n_discrete_pt -= sub->getNumberOfGhosts();
        return sub->getDomainID();
    }

    size_t getNumberOfSubDomains() const {return _list_dom.size();};
    DDCSubDomain* getSubDomain(size_t i) {return _list_dom[i];};

    size_t findDomainID(size_t global_obj_id) const
    {
        for (size_t i=0; i<_list_dom_start_obj.size(); ++i) {
            if (_list_dom_start_obj[i] <= global_obj_id)
                return i;
        }
        return 0; //error
    }

    size_t getTotalNumberOfDecomposedObjects() const {return _n_discrete_pt;}; 
private:
    DISALLOW_COPY_AND_ASSIGN(DDCGlobal);

private:
    DecompositionType::type _ddc_type;
    std::vector<size_t> _list_dom_start_obj;
    std::vector<DDCSubDomain*> _list_dom;
    size_t _n_discrete_pt;
    size_t _id;
};

} //end

