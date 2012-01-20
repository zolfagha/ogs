
#pragma once

namespace NumLib
{

template<typename Tvalue, typename Tpos>
class IFunction
{
public:
    virtual Tvalue& getValue(Tpos &pt) const = 0;
};

}
