
#pragma once

#include "BaseLib/CodingTools.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "IFeObjectContainer.h"
#include "FeObjectCachePerFeType.h"

namespace FemLib
{

/**
 * \brief Lagrangian finite element object containers
 */
class LagrangianFeObjectContainer : public FeObjectCachePerFeType, public IFeObjectContainer
{
public:
    LagrangianFeObjectContainer(MeshLib::IMesh &msh) :  FeObjectCachePerFeType(msh), _order(1) {};
    virtual ~LagrangianFeObjectContainer() {};

    void setPolynomialOrder(size_t order) { _order = order; };

    virtual IFiniteElement* getFeObject(const MeshLib::IElement &e);

private:
    FiniteElementType::type getFeType(MeshLib::ElementShape::type ele_type, size_t order);

    DISALLOW_COPY_AND_ASSIGN(LagrangianFeObjectContainer);

private:
    size_t _order;
};



}
