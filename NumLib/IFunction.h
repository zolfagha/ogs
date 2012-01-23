
#pragma once

#include "GeoLib/Point.h"

namespace NumLib
{

template<typename Tvalue>
class IFunction
{
public:
    virtual Tvalue& getValue(GeoLib::Point &pt) const = 0;
};

}
