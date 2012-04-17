
#pragma once

#include "IOSystem.h"
#include "NamedIO.h"
#include "MathLib/Parameter/DynamicIOSystem.h"

namespace MathLib
{

template <class T_BASE>
class NamedIOSystem : public NamedIO<IOSystem<T_BASE> >
{
public:
	virtual ~NamedIOSystem() {};
};

template <class T_BASE>
class NamedDynamicIOSystem : public NamedIO<DynamicIOSystem<T_BASE> >
{
public:
	virtual ~NamedDynamicIOSystem() {};
};

} //end
