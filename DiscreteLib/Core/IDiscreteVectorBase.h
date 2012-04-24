
#pragma once

#include <cstddef>

namespace DiscreteLib
{

/**
 * \brief Interface of Vector containers in discrete systems
 */
class IDiscreteVectorBase
{
public:
	IDiscreteVectorBase() : _obj_id(0) {};
    virtual ~IDiscreteVectorBase() {};
    size_t getObjectID() const {return _obj_id;};
    void setObjectID(size_t i) {_obj_id = i;};
private:
    size_t _obj_id;
};

} // end
