/*
 * GeoObject.h
 *
 *  Created on: Aug 27, 2010
 *      Author: TF
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

namespace GeoLib {

struct GeoObjType
{
    enum type {
        POINT,
        POLYLINE,
        POLYGON,
        SURFACE,
        VOLUME,
        INVALID
    };
};
    

/**
 * \ingroup GEOLIB
 *
 * \brief Base class for classes Point, Polyline, Surface.
 */

class GeoObject {
public:
    GeoObject() {};
    virtual ~GeoObject() {};

    virtual GeoObjType::type getGeoType() const = 0;
};

} // end namespace GEOLIB

#endif /* GEOOBJECT_H_ */
