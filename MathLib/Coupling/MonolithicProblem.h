
#pragma once

#include <cassert>
#include <vector>
#include <queue>
#include "Base/CodingTools.h"
#include "MathLib/Parameter/NamedIOSystem.h"
#include "ICoupledProblem.h"

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

class TemplateSteadyMonolithicSystem : public AbstractMonolithicSystem<ICoupledSystem>
{
};

//template <class T_SUPER, size_t N_IN, size_t N_OUT>
//class TemplateMonolithicSystem : public AbstractMonolithicSystem<T_SUPER>
//{
//public:
//    TemplateMonolithicSystem() {
//        // assume parameter id is assigned as in0, in1, ..., out0, out1
//        for (size_t i=0; i<N_IN; i++)
//            AbstractMonolithicSystem<T_SUPER>::registerInputParameter(i);
//        for (size_t i=0; i<N_OUT; i++)
//            AbstractMonolithicSystem<T_SUPER>::addOutputParameter(N_IN + i);
//    }
//};
//
//template <size_t N_IN, size_t N_OUT>
//class TemplateSteadyMonolithicSystem : public TemplateMonolithicSystem<ICoupledSystem, N_IN, N_OUT>
//{
//};

}
