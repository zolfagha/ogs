
#pragma once

#include <vector>

#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"
//#include "NumLib/Coupling/Algorithm/PartitionedAlgorithm.h"
#include "TransientCoupledSystem.h"
#include "TransientPartitionedAlgorithm.h"

namespace NumLib
{

/**
 * \brief Asynchronized partitioned problem
 */
class AsyncPartitionedSystem : public ITransientCoupledSystem
{
	typedef MathLib::IFunction Variable;
	typedef VariableContainer MyNamedVariableContainer;
	typedef ITransientPartitionedAlgorithm MyTransientPartitionedAlgorithm;
	typedef ITransientCoupledSystem MyTransientCoupledSystem;
	typedef ICoupledSystem MyCoupledSystem;
    typedef VariableMappingTable MyVariableMappingTable;
public:
    AsyncPartitionedSystem(MyTransientPartitionedAlgorithm &algo) : _algorithm(&algo)
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
    Variable* getParameter(size_t para_id) const { return _vars_t_n1.get(para_id); }
    template <class T>
    T* getParameter(size_t para_id) const
    {
        return static_cast<T*>(_vars_t_n1.get(para_id));
    }

    /// find subproblem
    int find(const MyTransientCoupledSystem &sub) const;

    /// check consistency
    bool check() const;


    /// add parameter without giving reference
    size_t addParameter(const std::string &name);

    /// add parameter and reference
    /// @param name variable name
    /// @param sys problem
    /// @param internal_id parameter id in the sys
    /// @return parameter id
    size_t addParameter(const std::string &name, MyTransientCoupledSystem & sub_problem, size_t para_id_in_sub_problem);

    /// set parameter 
    void setParameter(size_t para_id, Variable* var)
    {
        _vars_t_n1.set(para_id, *var);
    }

    void resetParameters( Variable* var)
    {
        for (size_t i=0; i<_vars_t_n1.size(); i++) {
            _vars_t_n1.set(i, *var);
        }
        _vars_t_n.clear();
        _vars_t_n1.clone(_vars_t_n);
    }

    /// connect system input and shared variable
    void connectInput(const std::string &this_para_name, MyTransientCoupledSystem &subproblem, size_t subproblem_para_id);

    //void addChildren(ITransientCoupledProblem& sys);
    //size_t getNumberOfChildren() const;
	
	double suggestNext(const TimeStep &time_current);
	
	int solveTimeStep(const TimeStep &time);
	
	bool isAwake(const TimeStep &time);

    void accept(const TimeStep &time);

    void getActiveProblems(const TimeStep &time, std::vector<MyCoupledSystem*> &list_active_problems);

private:
    MyTransientPartitionedAlgorithm *_algorithm;
    std::vector<MyTransientCoupledSystem*> _list_subproblems;
    std::vector<TimeStep> _list_synchronize_time;
	
    //std::vector<ITransientCoupledProblem*> _list_subproblems;
    std::vector<size_t> _list_input_parameters;
    MyNamedVariableContainer _vars_t_n;
    MyNamedVariableContainer _vars_t_n1;
    MyVariableMappingTable _map;

    size_t addSubProblem(MyTransientCoupledSystem &sub_problem);
};


}
