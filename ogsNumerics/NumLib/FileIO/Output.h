
#pragma once

#include "GeoLib/Core/GeoObject.h"

namespace NumLib
{


struct TimeType
{
enum type
{
    ALL_STEPS
};
};

class IOutput
{
public:
    IOutput(TimeType::type /*time_tpye*/) {};

    void addNodal(int);
    void addElemental(int);
};

class PVD : public IOutput
{
public:
    PVD(TimeType::type time_tpye) : IOutput(time_tpye) {};
};


class TECPLOT : public IOutput
{
public:
    TECPLOT(TimeType::type time_tpye) : IOutput(time_tpye) {};
    void setGeometry(GeoLib::GeoObject& geo);
};

}
