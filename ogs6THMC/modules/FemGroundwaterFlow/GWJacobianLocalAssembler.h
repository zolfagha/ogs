
#pragma once

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/LagrangeFeObjectContainer.h"
#include "NumLib/TransientAssembler/IElementWiseTransientJacobianLocalAssembler.h"
#include "MaterialLib/PorousMedia.h"


class GroundwaterFlowJacobianLocalAssembler
: public NumLib::IElementWiseTransientJacobianLocalAssembler
{
public:
    explicit GroundwaterFlowJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects)
    : _feObjects(&feObjects)
    {
    };

//    GroundwaterFlowJacobianLocalAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MaterialLib::PorousMedia &pm)
//    : _pm(&pm), _feObjects(&feObjects)
//    {
//    };

    void assembly(const NumLib::TimeStep &/*time*/, MeshLib::IElement &e, const NumLib::LocalVector &/*u1*/, const NumLib::LocalVector &/*u0*/,  NumLib::LocalMatrix &localJ);

private:
    //MaterialLib::PorousMedia* _pm;
    FemLib::LagrangianFeObjectContainer* _feObjects;
};


