
#pragma once

#include "MeshLib/Core/IElement.h"

#include "AbstractFemIntegrationGaussBase.h"
#include "FemIntegrationGaussLine.h"
#include "FemIntegrationGaussQuad.h"
#include "FemIntegrationGaussTriangle.h"

namespace FemLib
{

class FemGaussIntegrationFactory
{
public:
    static AbstractFemIntegrationGaussBase* create(MeshLib::IElement *e) 
    {
        switch (e->getShapeType()) {
            case MeshLib::ElementShape::QUAD:
                return new FemIntegrationGaussQuad();
            case MeshLib::ElementShape::TRIANGLE:
                return new FemIntegrationGaussQuad();
            default:
                return 0;
        }
    }
};

}
