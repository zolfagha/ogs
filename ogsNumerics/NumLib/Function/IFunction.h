
#pragma once

#include "IClonable.h"

namespace NumLib
{

/**
 *
 */
class IFunction : public IClonable
{
public:
    virtual ~IFunction() {};
};

/**
 *
 */
template<typename Tpos, typename Tval>
class TemplateFunction : public IFunction
{
public:
    virtual ~TemplateFunction() {};
    virtual void eval(const Tpos &x, Tval &val) = 0;
};

/**
 *
 */
template<typename Tpos, typename Tval, class T_FUNCTION>
class AbstractDefaultCloneFunction : public TemplateFunction<Tpos, Tval>
{
public:
    virtual ~AbstractDefaultCloneFunction() {};
    virtual IFunction* clone() const
    {
        return new T_FUNCTION();
    }
};

}
