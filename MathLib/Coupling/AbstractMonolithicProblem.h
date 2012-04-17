
#pragma once

#include <cassert>
#include <vector>
#include <queue>
#include "Base/CodingTools.h"
#include "MathLib/Parameter/NamedIOSystem.h"

namespace MathLib
{

/**
 * \brief MonolithicSolution
 */
template<class T_BASE>
class AbstractMonolithicSystem : public NamedIOSystem<T_BASE>
{
public:
	///
    AbstractMonolithicSystem() {};

    ///
    virtual ~AbstractMonolithicSystem()
    {
    }

    /// check consistency
    virtual bool check() const {return true;};

private:

    DISALLOW_COPY_AND_ASSIGN(AbstractMonolithicSystem);
};

}
