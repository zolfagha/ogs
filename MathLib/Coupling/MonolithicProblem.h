
#pragma once

#include <cassert>
#include <vector>

#include "Base/CodingTools.h"

#include "ICoupledProblem.h"

namespace MathLib
{

/**
 * \brief MonolithicSolution
 */
template<class T_SUPER>
class AbstractMonolithicSystem : public T_SUPER
{
public:
	typedef MathLib::IFunction Variable;

    virtual ~AbstractMonolithicSystem()
    {
        Base::releaseObjectsInStdVector(_vec_parameters);
    }

    void setParameter(size_t in_var, Variable* var)
    {
        assert(in_var<_vec_parameters.size());
        if (_vec_parameters[in_var]!=0)
            delete _vec_parameters[in_var];
        _vec_parameters[in_var] = var->clone();
    }

    Variable* getParameter(size_t out) const
    {
        assert(out<_vec_parameters.size());
        return _vec_parameters[out];
    }

    size_t getParameterIdForInput(size_t input_id) const {return input_id;};

    bool check() const {return true;};

protected:
    std::vector<MathLib::IFunction*> _vec_parameters;
//    std::vector<Variable*> _vec_parameters;
};

template <class T_SUPER, size_t N_IN, size_t N_OUT>
class TemplateMonolithicSystem : public AbstractMonolithicSystem<T_SUPER>
{
public:
    TemplateMonolithicSystem() {
	    AbstractMonolithicSystem<T_SUPER>::_vec_parameters.resize(getNumberOfParameters());
    }

    size_t getNumberOfInputParameters() const {return N_IN;};
    size_t getNumberOfParameters() const {return N_IN+N_OUT;};
};

template <size_t N_IN, size_t N_OUT>
class TemplateSteadyMonolithicSystem : public TemplateMonolithicSystem<ICoupledSystem, N_IN, N_OUT>
{
};

}
