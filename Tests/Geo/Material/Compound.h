/*
 * Compound.h
 *
 *  Created on: 10.04.2012
 *      Author: watanabe
 */

#ifndef COMPOUND_H_
#define COMPOUND_H_

#include "Base/CodingTools.h"

namespace NumLib
{
class ITXFunction;
}

namespace Geo
{

struct Compound
{
	NumLib::ITXFunction* molecular_diffusion;

	Compound()
	{
		Base::zeroObject(molecular_diffusion);
	}
	~Compound()
	{
		Base::releaseObject(molecular_diffusion);
	}

};

} //end

#endif /* COMPOUND_H_ */
