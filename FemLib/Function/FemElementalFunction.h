
#pragma once

#include <vector>
#include <algorithm>

#include "Base/CodingTools.h"

#include "MathLib/Vector.h"
#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/Core/DiscreteVector.h"

#include "NumLib/Function/Function.h"

#include "FemLib/Core/Element/IFemElement.h"
#include "FemLib/Tools/FemElementObjectContainer.h"

namespace FemLib
{



/**
 * \brief Template class for FEM element-based functions
 */
template<typename Tvalue>
class TemplateFEMElementalFunction : public NumLib::TemplateFunction<GeoLib::Point,Tvalue>
{
public:
    void eval(const GeoLib::Point &pt, Tvalue &v) {
        v = _ele_values[0];
    };
    Tvalue& getValue(size_t id) const {
        return _ele_values[id];
    };
private:
    Tvalue* _ele_values;
};


}



