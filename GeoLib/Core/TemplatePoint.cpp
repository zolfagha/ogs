
#include "TemplatePoint.h"

namespace GeoLib
{

#ifndef MSVC
template <>
bool TemplatePoint<double>::operator== (const TemplatePoint<double> &p) const {
    for (size_t i=0; i<3; i++)
        if (fabs(_x[i]-p._x[i])>std::numeric_limits<double>::epsilon()) return false;
    return true;
}

template class TemplatePoint<double>;
#endif

} //end

