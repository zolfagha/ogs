
#pragma once

#include <map>

#include "FemLib/IFemElement.h"
#include "FemLib/Element/LINE2.h"
#include "FemLib/Element/QUAD4.h"
#include "FemLib/Element/TRI3CONST.h"

namespace FemLib
{

class FemElementFactory
{
public:
    static IFiniteElement* create(FiniteElementType::type fe_type)
    {
        switch (fe_type)
        {
            case FiniteElementType::LINE2:
                return new LINE2();
            case FiniteElementType::QUAD4:
                return new QUAD4();
            case FiniteElementType::TRI3CONST:
                return new TRI3CONST();
        }
        assert(false);
        return 0;
    }
};

}
