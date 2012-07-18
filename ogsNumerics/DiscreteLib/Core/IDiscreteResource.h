
#pragma once

#include <cstddef>

namespace DiscreteLib
{

/**
 * \brief Interface of all resource
 */
class IDiscreteResource
{
public:
    IDiscreteResource() : _obj_id(0) {};
    virtual ~IDiscreteResource() {};
    size_t getObjectID() const {return _obj_id;};
    void setObjectID(size_t i) {_obj_id = i;};
private:
    size_t _obj_id;
};

} // end
