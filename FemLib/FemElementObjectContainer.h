
#pragma once

#include <map>
#include "Base/MemoryTools.h"
#include "FemLib/IFemElement.h"
#include "FemLib/Element/FemElementFactory.h"

namespace FemLib
{

class FeObjectCachePerFeType
{
public:
    virtual ~FeObjectCachePerFeType()
    {
        destroyStdMapWithPointers(_mapFeObj);
    }

    IFiniteElement* getFeObject(FiniteElementType::type fe_type)
    {
        IFiniteElement *fe = 0;
        if (_mapFeObj.count(fe_type)==0) {
            fe = FemElementFactory::create(fe_type);
            _mapFeObj[fe_type] = fe;
        } else {
            fe = _mapFeObj[fe_type];
        }
        return fe;
    }

private:
    std::map<FiniteElementType::type, IFiniteElement*> _mapFeObj;
};


class IFeObjectContainer
{
public:
    virtual IFiniteElement* getFeObject(MeshLib::IElement *e) = 0;
};

class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    FeObjectContainerPerElement(size_t ele_size)
    {
        _vec_fem.resize(ele_size);
    }
    virtual ~FeObjectContainerPerElement()
    {
        destroyStdVectorWithPointers(_vec_fem);
    }

    void addFiniteElement(size_t i, FiniteElementType::type fe_type)
    {
        _vec_fem[i] = FemElementFactory::create(fe_type);
    }

    virtual IFiniteElement* getFeObject(MeshLib::IElement *e) 
    {
        return _vec_fem[e->getID()];
    }

private:
    std::vector<IFiniteElement*> _vec_fem;
};



class LagrangianFeObjectContainer : public FeObjectCachePerFeType, public IFeObjectContainer
{
public:
    LagrangianFeObjectContainer() 
    {
        _order = 1;
    }

    void setPolynomialOrder(size_t order) 
    {
        _order = order;
    }

    virtual IFiniteElement* getFeObject(MeshLib::IElement *e)
    {
        FiniteElementType::type fe_type = getFeType(e->getElementType(), _order);
        return FeObjectCachePerFeType::getFeObject(fe_type);
    }
private:
    size_t _order;

    FiniteElementType::type getFeType(MeshLib::ElementType::type ele_type, size_t order)
    {
        switch (ele_type)
        {
            case MeshLib::ElementType::LINE:
                return (order==1) ? FiniteElementType::LINE2 : FiniteElementType::LINE3;
            case MeshLib::ElementType::QUAD:
                return (order==1) ? FiniteElementType::QUAD4 : FiniteElementType::QUAD9;
        }
        return FiniteElementType::INVALID;
    };
};



}
