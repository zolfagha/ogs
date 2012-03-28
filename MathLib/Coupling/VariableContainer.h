
#pragma once

#include <string>
//#include <map>
#include <vector>
#include <algorithm>

#include "Base/CodingTools.h"
#include "MathLib/Function/Function.h"

namespace MathLib
{

/**
 * \brief Container of named variables
 */
class VariableContainer
{
public:
	typedef MathLib::IFunction Variable;

	///
	VariableContainer() {};

    ///
    virtual ~VariableContainer()
    {
        Base::releaseObjectsInStdVector(_list_own_var_data);
    }

    /// reset data
    void clear()
    {
        _list_var_names.clear();
        _list_var_data.clear();
        Base::releaseObjectsInStdVector(_list_own_var_data);
    }

    /// make a copy of the given object
    void assign(const VariableContainer &src)
    {
    	this->clear();
        _list_var_names.assign(src._list_var_names.begin(), src._list_var_names.end());
        _list_var_data.resize(src._list_var_data.size());
        _list_own_var_data.resize(src._list_var_data.size());
        for (size_t i=0; i<src._list_var_data.size(); i++) {
            if (src._list_var_data[i]!=0)  {
                _list_var_data[i] = src._list_var_data[i]->clone();
                _list_own_var_data[i] = _list_var_data[i];
            }
        }
    }

    /// return the size of the stores variables
    size_t size() const 
    {
        return _list_var_names.size();
    }

    /// return if the container has a variable with the given key
    bool contain(const std::string &var_name) const
    {
    	return (this->find(var_name)>=0);
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

    /// return the variable with the given index
    template <class T>
    T* get(size_t var_id) const
    {
        return static_cast<T*>(_list_var_data[var_id]);
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
#if 0
    void set(size_t var_id,  const Variable& v)
    {
        _list_var_data[var_id] = v.clone();
    }
#else
    void set(size_t var_id,  Variable& v)
    {
        _list_var_data[var_id] = &v; //.clone();
    }
#endif

private:
    std::vector<std::string> _list_var_names;
    std::vector<Variable*> _list_var_data;
    std::vector<Variable*> _list_own_var_data;

    DISALLOW_COPY_AND_ASSIGN(VariableContainer);
};

}
