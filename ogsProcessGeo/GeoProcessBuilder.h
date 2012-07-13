
#pragma once

#include "ProcessLib/ProcessBuilder.h"

class GeoProcessBuilder: public ProcessLib::AbstractProcessBuilder
{
public:
	/// get a global instance
	static ProcessLib::AbstractProcessBuilder* getInstance();
private:
    static ProcessLib::AbstractProcessBuilder* _obj;

public:
	GeoProcessBuilder();
	virtual ~GeoProcessBuilder() {};
};
