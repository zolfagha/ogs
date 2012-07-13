
#include "GeoProcessBuilder.h"

#include "ProcessLib/ProcessFactoryImpl.h"
#include "ProcessList.h"

ProcessLib::AbstractProcessBuilder* GeoProcessBuilder::_obj = 0;

ProcessLib::AbstractProcessBuilder* GeoProcessBuilder::getInstance()
{
	if (_obj==0) _obj = new GeoProcessBuilder();
	return _obj;
}


GeoProcessBuilder::GeoProcessBuilder()
{
#include "ProcessReg.txt"
}
