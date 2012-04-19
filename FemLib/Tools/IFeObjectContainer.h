
#pragma once

#include "FemLib/Core/Element/IFemElement.h"

namespace FemLib
{

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

}
