
#pragma once

#include "ProcessLib/AbstractProcessBuilder.h"

/**
 * \brief Process builder for OGS6
 *
 * This class follows singleton pattern.
 */
class GeoProcessBuilder: public ProcessLib::AbstractProcessBuilder
{
public:
	static ProcessLib::AbstractProcessBuilder* getInstance();
private:
    static ProcessLib::AbstractProcessBuilder* _obj;

public:
	virtual ~GeoProcessBuilder() {};

private:
	GeoProcessBuilder();
};
