

#pragma once

#include "MeshLib/Core/IElement.h"
#include "FemLib/IFemElement.h"
#include "FemLib/FemElementObjectContainer.h"

namespace FemLib
{

class FemLagrangeElement : public IFiniteElement
{
private:

public:
    /// initialize object for given mesh elements
    virtual void configure( const MeshLib::IElement * e )
    {

    }
};

}
