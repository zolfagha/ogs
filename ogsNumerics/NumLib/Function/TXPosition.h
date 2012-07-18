
#pragma once

#include "GeoLib/Core/Point.h"

namespace NumLib
{

/**
 * \brief Position in space-time domain
 */
class TXPosition
{
public:
    TXPosition(double t_, const double* x_) : _t(t_), _x(x_), _id(0) {};
    TXPosition(const double* x_) : _t(.0), _x(x_), _id(0) {};
    TXPosition(const GeoLib::Point* p) : _t(.0), _x(p->getData()), _id(0) {};
    TXPosition(size_t n, const GeoLib::Point* p) : _t(.0), _x(p->getData()), _id(n) {};
    TXPosition(size_t n) : _t(.0), _x(0), _id(n) {};
    

    double getTime() const {return _t;};
    const double* getSpace() const {return _x;};
    size_t getId() const {return _id;};
private:
    double _t;
    const double* _x;
    size_t _id;
};


} //end
