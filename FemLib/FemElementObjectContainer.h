
#pragma once

#include <map>

#include "Base/CodingTools.h"

#include "FemLib/Core/IFemElement.h"
#include "FemLib/Core/Element/FemElementFactory.h"


namespace FemLib
{

/**
 * \brief Cache system for finite element objects 
 */
class FeObjectCachePerFeType
{
public:
    FeObjectCachePerFeType() {};
    virtual ~FeObjectCachePerFeType()
    {
        Base::releaseObjectsInStdMap(_mapFeObj);
    }

    IFiniteElement* getFeObject(FiniteElementType::type fe_type, MeshLib::IMesh &msh)
    {
        IFiniteElement *fe = 0;
        if (_mapFeObj.count(fe_type)==0) {
            fe = FemElementFactory::create(fe_type, msh);
            _mapFeObj[fe_type] = fe;
        } else {
            fe = _mapFeObj[fe_type];
        }
        return fe;
    }

private:
    std::map<FiniteElementType::type, IFiniteElement*> _mapFeObj;

    DISALLOW_COPY_AND_ASSIGN(FeObjectCachePerFeType);
};

/**
 * \brief Interface of finite element object containers
 */
class IFeObjectContainer
{
public:
    /// get a finite element object for the given mesh element
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, MeshLib::IMesh &msh) = 0;
};

/**
 * \brief Finite element object containers
 */
class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    FeObjectContainerPerElement(size_t ele_size)
    {
        _vec_fem.resize(ele_size);
    }
    virtual ~FeObjectContainerPerElement()
    {
        Base::releaseObjectsInStdVector(_vec_fem);
    }

    void addFiniteElement(size_t i, FiniteElementType::type fe_type, MeshLib::IMesh &msh)
    {
        _vec_fem[i] = FemElementFactory::create(fe_type, msh);
    }

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, MeshLib::IMesh &msh) 
    {
        return _vec_fem[e.getID()];
    }

private:
    std::vector<IFiniteElement*> _vec_fem;
    DISALLOW_COPY_AND_ASSIGN(FeObjectContainerPerElement);
};


/**
 * \brief Lagrangian finite element object containers
 */
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

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e, MeshLib::IMesh &msh)
    {
        FiniteElementType::type fe_type = getFeType(e.getShapeType(), _order);
        return FeObjectCachePerFeType::getFeObject(fe_type, msh);
    }

private:
    size_t _order;

    FiniteElementType::type getFeType(MeshLib::ElementShape::type ele_type, size_t order)
    {
        switch (ele_type)
        {
            case MeshLib::ElementShape::LINE:
                return (order==1) ? FiniteElementType::LINE2 : FiniteElementType::LINE3;
            case MeshLib::ElementShape::QUAD:
                return (order==1) ? FiniteElementType::QUAD4 : FiniteElementType::QUAD9;
        }
        return FiniteElementType::INVALID;
    };
    DISALLOW_COPY_AND_ASSIGN(LagrangianFeObjectContainer);
};



}
