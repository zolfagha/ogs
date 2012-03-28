
#pragma once

#include <vector>
#include <iostream>

#include "ICoupledProblem.h"
#include "ParameterTable.h"
#include "ParameterProblemMappingTable.h"

#include "MathLib/Coupling/Algorithm/PartitionedAlgorithm.h"

namespace MathLib
{

/**
 * \brief Partitioned problem
 */
class PartitionedProblem : public ICoupledSystem
{
public:
    /// 
    PartitionedProblem(IPartitionedAlgorithm &algo) : _algorithm(&algo)
    {
    }

    virtual ~PartitionedProblem() {};

    /// get the number of parameters
    size_t getNumberOfParameters() const {return _parameter_table.size();};

    /// get the number of input parameters
    size_t getNumberOfInputParameters() const {return _list_input_parameters.size();};

    /// get the index of parameter with the given name
    int getParameterID(const std::string &name) const { return _parameter_table.find(name); }

    /// get a id list of input parameters
    size_t getParameterIdForInput(size_t input_id) const {return _list_input_parameters[input_id];};

    /// get output parameter with the given id
    const Parameter* getOutput(size_t para_id) const { return _parameter_table.get(para_id); }

    /// get output parameter with the given id
    template <class T>
    const T* getOutput(size_t para_id) const
    {
    	return static_cast<const T*>(_parameter_table.get(para_id));
    }

    /// set output parameter
    void setOutput(size_t parameter_id, Parameter* val)
    {
        _parameter_table.set(parameter_id, *val);
    }

    /// set input parameter
    void setInput(size_t parameter_id, const Parameter* val)
    {
        _parameter_table.set(parameter_id, *val);
    }

    /// find subproblem
    int find(const ICoupledSystem& sub) const;

    /// check consistency
    bool check() const;


    /// add parameter without giving reference
    size_t addParameter(const std::string &name);

    /// add parameter and reference
    /// @param name variable name
    /// @param sys problem
    /// @param internal_id parameter id in the sys
    /// @return parameter id
    size_t addParameter(const std::string &name, ICoupledSystem &sub_problem, size_t para_id_in_sub_problem);

    /// connect system input and shared variable
    void connectInput(const std::string &this_para_name, ICoupledSystem &subproblem, size_t subproblem_para_id);

    /// solve this system
    int solve();


protected:
    /// get input parameter
    const Parameter* getInput(size_t parameter_id) const
    {
        return _parameter_table.get(parameter_id);
    }

private:
    std::vector<ICoupledSystem*> _list_subproblems;
    std::vector<size_t> _list_input_parameters;
    ParameterTable _parameter_table;
    ParameterProblemMappingTable _map;
    IPartitionedAlgorithm *_algorithm;

    size_t addSubProblem(ICoupledSystem &sub_problem);

};



}
