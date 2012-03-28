
#pragma once

#include <vector>

#include "MathLib/Coupling/ICoupledProblem.h"
#include "MathLib/Coupling/MonolithicProblem.h"
#include "MathLib/Coupling/PartitionedProblem.h"
#include "MathLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
//#include "NumLib/Coupling/Algorithm/PartitionedAlgorithm.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

/**
 * \brief Asynchronized partitioned problem
 */
class AsyncPartitionedSystem : public ITransientCoupledSystem
{
public:
    AsyncPartitionedSystem(MathLib::ITransientPartitionedAlgorithm &algo) : _algorithm(&algo)
    {
    }

    virtual ~AsyncPartitionedSystem() {};

    /// get the number of parameters
    size_t getNumberOfParameters() const {return _vars_t_n1.size();};

    /// get the number of input parameters
    size_t getNumberOfInputParameters() const {return _list_input_parameters.size();};

    /// get the index of parameter with the given name
    int getParameterID(const std::string &name) const { return _vars_t_n1.find(name); }

    size_t getParameterIdForInput(size_t input_id) const {return _list_input_parameters[input_id];};

    /// get parameter with the given id
    const MathLib::Parameter* getOutput(size_t para_id) const { return _vars_t_n1.get(para_id); }
    template <class T>
    const T* getOutput(size_t para_id) const
    {
        return static_cast<const T*>(_vars_t_n1.get(para_id));
    }

    /// find subproblem
    int find(const ITransientCoupledSystem &sub) const;

    /// check consistency
    bool check() const;


    /// add parameter without giving reference
    size_t addParameter(const std::string &name);

    /// add parameter and reference
    /// @param name variable name
    /// @param sys problem
    /// @param internal_id parameter id in the sys
    /// @return parameter id
    size_t addParameter(const std::string &name, ITransientCoupledSystem & sub_problem, size_t para_id_in_sub_problem);

    /// set parameter 
    //void setParameter(size_t para_id, MathLib::Parameter* var)
    //{
    //    _vars_t_n1.set(para_id, *var);
    //}

    void resetParameters( MathLib::Parameter* var)
    {
        for (size_t i=0; i<_vars_t_n1.size(); i++) {
            _vars_t_n1.set(i, *var);
        }
        _vars_t_n.assign(_vars_t_n1);
    }

    /// connect system input and shared variable
    void connectInput(const std::string &this_para_name, ITransientCoupledSystem &subproblem, size_t subproblem_para_id);

    //void addChildren(ITransientCoupledProblem& sys);
    //size_t getNumberOfChildren() const;
	
	double suggestNext(const TimeStep &time_current);
	
	int solveTimeStep(const TimeStep &time);
	
	bool isAwake(const TimeStep &time);

    void accept(const TimeStep &time);

    void getActiveProblems(const TimeStep &time, std::vector<ICoupledSystem*> &list_active_problems);

    void setOutput(size_t para_id, MathLib::Parameter* var)
    {
        _vars_t_n1.set(para_id, *var);
    }
    void setInput(size_t para_id, const MathLib::Parameter* var)
    {
        _vars_t_n1.set(para_id, *var);
    }
protected:
    const MathLib::Parameter* getInput(size_t parameter_id) const
    {
        return _vars_t_n1.get(parameter_id);
    }

private:
    MathLib::ITransientPartitionedAlgorithm *_algorithm;
    std::vector<ITransientCoupledSystem*> _list_subproblems;
    std::vector<TimeStep> _list_synchronize_time;
	
    //std::vector<ITransientCoupledProblem*> _list_subproblems;
    std::vector<size_t> _list_input_parameters;
    MathLib::ParameterTable _vars_t_n;
    MathLib::ParameterTable _vars_t_n1;
    MathLib::ParameterProblemMappingTable _map;

    size_t addSubProblem(ITransientCoupledSystem &sub_problem);
};


}
