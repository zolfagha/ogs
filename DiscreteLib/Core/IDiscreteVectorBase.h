
#pragma once

#include <cstddef>
#include "IDiscreteResource.h"

namespace DiscreteLib
{

/**
 * \brief Interface of Vector containers in discrete systems
 */
class IDiscreteVectorBase : public IDiscreteResource
{
public:
	IDiscreteVectorBase() {};
    virtual ~IDiscreteVectorBase() {};
};

} // end
