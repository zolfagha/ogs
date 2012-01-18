
#pragma once

#include "Point.h"

namespace GeoLib
{

/**
 * Checks if two points are within a given distance of each other
 * @param p0 The first point
 * @param p1 the second point
 * @param squaredDistance The square of the distance within which the two points should be
 * @return true if p1 and p2 are within the given distance of each other, false otherwise
 */
bool checkDistance(GeoLib::Point const &p0, GeoLib::Point const &p1, double squaredDistance);

/** squared euklid norm of the vector p0 */
double sqrNrm2(const GeoLib::Point* const p0);

/** squared dist between GEOLIB::Points p0 and p1 */
double sqrDist(const GeoLib::Point* p0, const GeoLib::Point* p1);


}