
#pragma once

#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "Base/MemoryTools.h"
#include "MathLib/Function/Function.h"

namespace NumLib
{

typedef MathLib::IFunction<double, double> Variable;  //TODO general function?


/**
 * \brief Container of named variables
 */
class NamedVariableContainer
{
public:
    /// destructor
    virtual ~NamedVariableContainer()
    {
        Base::destroyStdVectorWithPointers(_list_var_data);
    }

    /// reset data
    void clear()
    {
        _list_var_names.clear();
        Base::destroyStdVectorWithPointers(_list_var_data);
    }

    /// make a copy of this object
    /// @param dest the destination object
    void clone(NamedVariableContainer &dest) const
    {
        dest.clear();
        dest._list_var_names.assign(_list_var_names.begin(), _list_var_names.end());
        dest._list_var_data.resize(_list_var_data.size());
        for (size_t i=0; i<_list_var_data.size(); i++) {
            if (_list_var_data[i]!=0)  {
                dest._list_var_data[i] = _list_var_data[i]->clone();
            }
        }
    }

    /// return the size of the stores variables
    size_t size() const 
    {
        return _list_var_names.size();
    }

    /// return if the container has a variable with the given key
    size_t contain(const std::string &var_name) const
    {
        if (find(var_name)<0)
            return 0;
        else
            return 1;
    }

    /// return the name of the variable with the given index
    const std::string& getName(size_t var_id) const 
    {
        return _list_var_names[var_id];
    };

    /// return the index of the variable with the given name
    int find(const std::string& name) const 
    {
        std::vector<std::string>::const_iterator itr = std::find(_list_var_names.begin(), _list_var_names.end(), name);
        if (itr!=_list_var_names.end()) {
            return (itr - _list_var_names.begin());
        }

        return -1;
    }

    /// return the variable with the given index
    Variable* get(size_t var_id) const
    {
        return _list_var_data[var_id];
    }

    /// register the variable
    size_t add(const std::string &var_name)
    {
        size_t new_id = _list_var_names.size();
        _list_var_names.push_back(var_name);
        _list_var_data.push_back(0);
        return new_id;
    }

    /// set variable
    void set(size_t var_id,  const Variable& v)
    {
        _list_var_data[var_id] = v.clone();
    }

private:
    std::vector<std::string> _list_var_names;
    std::vector<Variable*> _list_var_data;
};

class ICoupledProblem;

class VariableMappingTable
{
public:
    typedef std::pair<ICoupledProblem*,size_t> PairSysVarId;
    typedef std::pair<size_t, size_t> PairInputVar;
    typedef std::vector<PairInputVar> ListOfInputVar;
    std::vector<ListOfInputVar> _list_subproblem_input_source;
    std::vector<PairSysVarId> _map_paraId2subproblem;
};


}
