
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
	TXPosition(double t_, const double* x_) : _t(t_), _x(x_) {};
	TXPosition(const double* x_) : _x(x_) {};
	TXPosition(const GeoLib::Point* p) : _x(p->getData()) {};

	double getTime() const {return _t;};
	const double* getSpace() const {return _x;};
private:
	double _t;
	const double* _x;
};


} //end
