
#pragma once

#include <map>

#include "Base/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
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

    IFiniteElement* getFeObject(FiniteElementType::type fe_type);

    MeshLib::IMesh* getMesh() const {return _msh;};

private:
    std::map<FiniteElementType::type, IFiniteElement*> _mapFeObj;
    MeshLib::IMesh *_msh;

    DISALLOW_COPY_AND_ASSIGN(FeObjectCachePerFeType);
};

}
