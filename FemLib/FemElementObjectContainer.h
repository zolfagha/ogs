
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
    FeObjectCachePerFeType(MeshLib::IMesh &msh) : _msh(&msh) {};

    virtual ~FeObjectCachePerFeType()
    {
        Base::releaseObjectsInStdMap(_mapFeObj);
    }

    IFiniteElement* getFeObject(FiniteElementType::type fe_type)
    {
        IFiniteElement *fe = 0;
        if (_mapFeObj.count(fe_type)==0) {
            fe = FemElementFactory::create(fe_type, *_msh);
            _mapFeObj[fe_type] = fe;
        } else {
            fe = _mapFeObj[fe_type];
        }
        return fe;
    }

    MeshLib::IMesh* getMesh() const {return _msh;};

private:
    std::map<FiniteElementType::type, IFiniteElement*> _mapFeObj;
    MeshLib::IMesh *_msh;

    DISALLOW_COPY_AND_ASSIGN(FeObjectCachePerFeType);
};

/**
 * \brief Interface of finite element object containers
 */
class IFeObjectContainer
{
public:
	virtual ~IFeObjectContainer() {};
    /// get a finite element object for the given mesh element
    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e) = 0;
};

/**
 * \brief Finite element object containers
 */
class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    FeObjectContainerPerElement(MeshLib::IMesh &msh, size_t ele_size) : _msh(&msh)
    {
        _vec_fem.resize(ele_size);
    }
    virtual ~FeObjectContainerPerElement()
    {
        Base::releaseObjectsInStdVector(_vec_fem);
    }

    void addFiniteElement(size_t i, FiniteElementType::type fe_type)
    {
        _vec_fem[i] = FemElementFactory::create(fe_type, *_msh);
    }

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e) 
    {
        return _vec_fem[e.getID()];
    }

private:
    MeshLib::IMesh* _msh;
    std::vector<IFiniteElement*> _vec_fem;
    DISALLOW_COPY_AND_ASSIGN(FeObjectContainerPerElement);
};


/**
 * \brief Lagrangian finite element object containers
 */
class LagrangianFeObjectContainer : public FeObjectCachePerFeType, public IFeObjectContainer
{
public:
    LagrangianFeObjectContainer(MeshLib::IMesh &msh) :  FeObjectCachePerFeType(msh)
    {
        _order = 1;
    }

    void setPolynomialOrder(size_t order) 
    {
        _order = order;
    }

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e)
    {
        FiniteElementType::type fe_type = getFeType(e.getShapeType(), _order);
        IFiniteElement* fe = FeObjectCachePerFeType::getFeObject(fe_type);
        fe->configure(*const_cast<MeshLib::IElement*>(&e));
        return fe;
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
            default:
                return FiniteElementType::INVALID;
        }
    };
    DISALLOW_COPY_AND_ASSIGN(LagrangianFeObjectContainer);
};



}
