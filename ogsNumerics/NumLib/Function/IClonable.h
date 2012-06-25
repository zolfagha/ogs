
#pragma once

namespace NumLib
{

class IClonable
{
public:
    virtual ~IClonable() {};
    virtual IClonable* clone() const = 0;
};

}
