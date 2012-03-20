#ifndef CPUTIMETIMER_H
#define CPUTIMETIMER_H

#include <ctime>

#include "TimeMeasurementBase.h"

namespace Base
{

class CPUTimeTimer
{
public:
    virtual void start();
    virtual void stop();
    virtual double elapsed();
	virtual ~CPUTimeTimer() {};
private:
	clock_t _start;
	clock_t _stop;
};

}

#endif

