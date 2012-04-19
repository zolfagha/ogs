
#pragma once

namespace MathLib
{

class IClonable
{
public:
    virtual ~IClonable() {};
    virtual IClonable* clone() const = 0;
};

}
