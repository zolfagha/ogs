
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
    virtual ~AbstractMonolithicSystem()
    {
    }

    void setNumberOfParameters(size_t n)
    {
        _vec_in_parameters.resize(n, 0);
        _vec_out_parameters.resize(n, 0);
    }

    void setInput(size_t parameter_id, const Parameter* val)
    {
        assert(parameter_id<_vec_in_parameters.size());
        _vec_in_parameters[parameter_id] = val;
    }

    const Parameter* getOutput(size_t parameter_id) const
    {
        assert(parameter_id<_vec_out_parameters.size());
        return _vec_out_parameters[parameter_id];
    }

    template <class T>
    const T* getOutput(size_t parameter_id) const
    {
        assert(parameter_id<_vec_out_parameters.size());
        return static_cast<T*>(_vec_out_parameters[parameter_id]);
    }

    size_t getParameterIdForInput(size_t parameter_id) const {return parameter_id;};

    bool check() const {return true;};

    void setOutput(size_t parameter_id, Parameter* val)
    {
        assert(parameter_id<_vec_out_parameters.size());
        _vec_out_parameters[parameter_id] = val;
    }
protected:
    const Parameter* getInput(size_t parameter_id) const
    {
        return _vec_in_parameters[parameter_id];
    }

    template <class T>
    const T* getInput(size_t parameter_id) const
    {
        assert(parameter_id<_vec_in_parameters.size());
        return static_cast<const T*>(_vec_in_parameters[parameter_id]);
    }

private:
    std::vector<const Parameter*> _vec_in_parameters;
    std::vector<Parameter*> _vec_out_parameters;
};

template <class T_SUPER, size_t N_IN, size_t N_OUT>
class TemplateMonolithicSystem : public AbstractMonolithicSystem<T_SUPER>
{
public:
    TemplateMonolithicSystem() {
	    AbstractMonolithicSystem<T_SUPER>::setNumberOfParameters(getNumberOfParameters());
    }

    size_t getNumberOfInputParameters() const {return N_IN;};
    size_t getNumberOfParameters() const {return N_IN+N_OUT;};
};

template <size_t N_IN, size_t N_OUT>
class TemplateSteadyMonolithicSystem : public TemplateMonolithicSystem<ICoupledSystem, N_IN, N_OUT>
{
};

}
