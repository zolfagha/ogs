/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplatePoint.cpp
 *
 * Created on 2010-01-28 by Thomas Fischer
 */

#include "TemplatePoint.h"

namespace GeoLib
{

#if 0
template <>
bool TemplatePoint<double>::operator== (const TemplatePoint<double> &p) const {
    for (size_t i=0; i<3; i++)
        if (fabs(_x[i]-p._x[i])>std::numeric_limits<double>::epsilon()) return false;
    return true;
}

template class TemplatePoint<double>;
#endif

} //end

