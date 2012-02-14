
#pragma once

#include <vector>

#include "ICoupledProblem.h"
#include "NamedVariableContainer.h"


namespace NumLib
{

class IPartitionedAlgorithm;

/**
 * \brief Partitioned problem
 */
class PartitionedProblem : public ICoupledProblem
{
public:

    /// 
    PartitionedProblem(IPartitionedAlgorithm &algo) : _algorithm(&algo)
    {
    }

    /// get the number of parameters
    size_t getNumberOfParameters() const {return _vars.size();};

    /// get the number of input parameters
    size_t getNumberOfInputParameters() const {return _list_input_parameters.size();};

    /// get the index of parameter with the given name
    int getParameterID(const std::string &name) const { return _vars.find(name); }

    size_t getParameterIdForInput(size_t input_id) const {return _list_input_parameters[input_id];};

    /// get parameter with the given id
    Variable* getParameter(size_t para_id) const { return _vars.get(para_id); }

    /// find subproblem
    int find(const ICoupledProblem& sub) const;

    /// check consistency
    bool check() const;


    /// add parameter without giving reference
    size_t addParameter(const std::string &name);

    /// add parameter and reference
    /// @param name variable name
    /// @param sys problem
    /// @param internal_id parameter id in the sys
    /// @return parameter id
    size_t addParameter(const std::string &name, ICoupledProblem& sub_problem, size_t para_id_in_sub_problem);

    /// set parameter 
    void setParameter(size_t para_id, Variable* var)
    {
        _vars.set(para_id, *var);
    }

    /// connect system input and shared variable
    void connectInput(const std::string &this_para_name, ICoupledProblem &subproblem, size_t subproblem_para_id);

    /// solve this system
    int solve();

private:
    std::vector<ICoupledProblem*> _list_subproblems;
    std::vector<size_t> _list_input_parameters;
    NamedVariableContainer _vars;
    VariableMappingTable _map;
    IPartitionedAlgorithm *_algorithm;

    size_t addSubProblem(ICoupledProblem &sub_problem);

};


}