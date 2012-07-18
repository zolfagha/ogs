/*
 * Compound.h
 *
 *  Created on: 10.04.2012
 *      Author: watanabe
 */

#ifndef COMPOUND_H_
#define COMPOUND_H_

#include "BaseLib/CodingTools.h"

namespace NumLib
{
class ITXFunction;
}

namespace MaterialLib
{

struct Compound
{
    NumLib::ITXFunction* molecular_diffusion;

    Compound()
    {
        BaseLib::zeroObject(molecular_diffusion);
    }
    ~Compound()
    {
        BaseLib::releaseObject(molecular_diffusion);
    }

};

} //end

#endif /* COMPOUND_H_ */
