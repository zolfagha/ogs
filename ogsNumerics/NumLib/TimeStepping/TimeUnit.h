
#pragma once

namespace NumLib
{

struct TimeUnit
{
    enum type
    {
        Second,
        Minute,
        Hour,
        Day,
        Week,
        Year
    };

    static double convertToSec(TimeUnit::type type, double v)
    {
        double fac = 1.0;
        switch (type)
        {
        case TimeUnit::Minute:
            fac = 60.0;
            break;
        case TimeUnit::Hour:
            fac = 3600.0;
            break;
        case TimeUnit::Day:
            fac = 86400.;
            break;
        case TimeUnit::Week:
            fac = 86400.*7.0;
            break;
        case TimeUnit::Year:
            fac = 86400.*365.0;
            break;
        default:
            break;
        }
        return v*fac;
    }
};

}
