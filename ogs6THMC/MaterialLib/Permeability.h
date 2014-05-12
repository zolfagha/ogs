/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Function/TXPosition.h"
#include "logog.hpp"


namespace MaterialLib
{

class IPorosityHydraulicConductivityModel
{
public:
	virtual void eval(double porosity, double &K) const = 0;
	virtual void eval(MathLib::LocalMatrix &porosity, MathLib::LocalMatrix &K) const
	{
		double k = .0;
		eval(porosity(0,0), k);
		K.resize(1,1);
		K(0,0) = k;
	}
	virtual ~IPorosityHydraulicConductivityModel() {}
	virtual IPorosityHydraulicConductivityModel* clone() const = 0;
};

class PorosityHydraulicConductivityModel1 : public IPorosityHydraulicConductivityModel
{
public:
	explicit PorosityHydraulicConductivityModel1(double porosity0, double K0) : _porosity0(porosity0), _K0(K0) {};

	virtual void eval(double porosity, double &K) const OGS_DECL_OVERRIDE
	{
		K = porosity/_porosity0*_K0;
	}

	virtual ~PorosityHydraulicConductivityModel1() {}

	virtual PorosityHydraulicConductivityModel1* clone() const OGS_DECL_OVERRIDE
	{
		return new PorosityHydraulicConductivityModel1(_porosity0, _K0);
	}

private:
	const double _porosity0;
	const double _K0;
};


} //end namespace
 
