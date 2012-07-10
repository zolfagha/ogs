
#include "ProcessBuilder.h"

#ifndef PROCRSS_REGISTER
#include "ProcessFactoryImpl.h"
#include "ProcessLib/ProcessList.h"

namespace ProcessLib
{
ProcessBuilder::ProcessBuilder()
{
//	this->registerProcess("GW", new ProcessFactoryImpl<FunctionHead>);
#include "ProcessLib/ProcessReg.txt"
}
}

#endif
