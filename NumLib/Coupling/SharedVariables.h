
#pragma once

#include <string>
#include <map>
#include <vector>

#include "MathLib/Function/Function.h"

namespace NumLib
{

typedef MathLib::IFunction<double, double> Variable;  //TODO general function?

class SharedVariables
{
private:
    std::vector<std::string> _list_var_names;
    std::map<std::string, Variable*> _map_id_func;

public:
    size_t registerVarialbe(const std::string &var_name)
    {
        size_t new_id = _list_var_names.size();
        //_map_id_func[var_name] = var;
        _list_var_names.push_back(var_name);
        return new_id;
    }

    size_t size() const 
    {
        return _list_var_names.size();
    }

    bool hasVariable(const std::string &var_name) const
    {
        return (_list_var_names.end()!=std::find(_list_var_names.begin(), _list_var_names.end(), var_name));
    }

    const std::string& getVariableName(size_t i) const 
    {
        return _list_var_names[i];
    };

    int getVariableID(const std::string& name) const 
    {
        for (size_t i=0; i<_list_var_names.size(); i++) {
            if (_list_var_names[i].compare(name)==0) return i;
        }
        return -1;
    }

    void setVariable(const std::string &var_name, Variable* v)
    {
        _map_id_func[var_name] = v;
    }

    void setVariable(size_t var_id,  Variable* v)
    {
        setVariable(getVariableName(var_id), v);
    }

    Variable* getVariable(const std::string &var_name) const
    {
        std::map<std::string, Variable*>::const_iterator itr = _map_id_func.find(var_name);
        if (itr==_map_id_func.end())
            return 0;
        else
            return itr->second;
    }

    Variable* getVariable(size_t var_id) const
    {
        return getVariable(getVariableName(var_id));
    }

    void reset()
    {
        _map_id_func.clear();
    }
};

}
