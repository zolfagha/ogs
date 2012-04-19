
#pragma once

#include "MathLib/MathTools.h"

inline double analyticalOgataBank(double x, double t, double v, double alpha)
{
    double u = MathLib::erfc((x-v*t)/std::sqrt(4*alpha*t))+std::exp(v*x/alpha)*MathLib::erfc((x+v*t)/std::sqrt(4*alpha*t));
    u *= 0.5;
    return u;
}

