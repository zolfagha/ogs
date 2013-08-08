/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeStepFunctionNewtonAdaptive.cpp
 *
 * Created on 2013-08-07 by Haibing Shao
 */

#include "TimeStepFunctionNewtonAdaptive.h"


namespace NumLib
{


TimeStepFunctionNewtonAdaptive::TimeStepFunctionNewtonAdaptive(double t0,
							   double tn,
							   double min_ts,
							   double max_ts,
							   const std::vector<int> &iter_times_vector,
							   const std::vector<double> &multiplier_vector,
							   TimeUnit::type t_unit)
: AbstractTimeStepFunction(t0, tn, t_unit), _min_ts( min_ts ), _max_ts( max_ts ), _iter_times(0)
{
	int tmp_iter_times;
	double tmp_multiplier;
	for (size_t i=0; i < iter_times_vector.size(); i++ )
	{
		tmp_iter_times = iter_times_vector[i];
		tmp_multiplier = multiplier_vector[i];
		_iter_times_vector.push_back( tmp_iter_times );
		_multiplier_vector.push_back( tmp_multiplier );
	}
};

TimeStepFunctionNewtonAdaptive* TimeStepFunctionNewtonAdaptive::clone()
{
	TimeStepFunctionNewtonAdaptive *obj = new TimeStepFunctionNewtonAdaptive(getBeginning(),
																			 getEnd(),
																			 _min_ts,
																			 _max_ts,
																			 _iter_times_vector,
																			 _multiplier_vector);
	return obj;
}

void TimeStepFunctionNewtonAdaptive::updateLog(BaseLib::Options &log)
{
	BaseLib::Options* optNR = log.getSubGroup("NewtonRaphson");
   if (optNR==nullptr) return;

   this->_iter_times = optNR->getOptionAsNum<size_t>("iterations");
};

double TimeStepFunctionNewtonAdaptive::suggestNext(double /*t_current*/)
{
	double t=0.0;

	if ( getStep() == 0 )  // the first time step
	{
		t = getPrevious() + _min_ts;
	}
	else  // not the first time step
	{
		double tmp_dt = 0.0;
		double tmp_multiplier = 1.0;
		double t_pre  = getPrevious();

		// number of iterations must be provided externally from the solver class
		int iter_times = this->_iter_times;  // to be changes.
		// finding the right multiplier
		for ( size_t i=0; i<_iter_times_vector.size(); i++ )
			if ( iter_times > _iter_times_vector[i] )
				tmp_multiplier = _multiplier_vector[i];
		// times the muliplier
		tmp_dt = t_pre * tmp_multiplier;
		// check whether out of the boundary
		if ( tmp_dt < _min_ts )
			tmp_dt = _min_ts;
		else if ( tmp_dt > _max_ts )
			tmp_dt = _max_ts;
		// getting the next time step point
		t = getPrevious() + tmp_dt;
	}

	if (t > getEnd())
		return getEnd();
	else
		return t;
};


}
