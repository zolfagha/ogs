
#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include "CouplingSolution.h"
#include "PartitionedAlgorithm.h"

namespace NumLib
{

/**
 * \brief Partitioned solution
 */
class PartitionedProblem : public ICouplingProblem
{
private:
    IPartitionedAlgorithm *_algorithm;
    std::vector<ICouplingProblem*> _listChildren;
    SharedVariables _vars;
    typedef std::pair<ICouplingProblem*,size_t> PairSysVarId;
    std::vector<std::pair<std::string, PairSysVarId>> _map_sharedVar2sysVar;
    typedef std::pair<size_t, std::string> PairInputVar;
    std::map<ICouplingProblem*, std::vector<PairInputVar>> _map_sys2input;
    std::vector<size_t> _list_input_var;

public:
    PartitionedProblem(IPartitionedAlgorithm *algo)
    {
        _algorithm = algo;
    }

    /// add shared variable and associated system
    void add(const std::string &var, ICouplingProblem* sys, int internal_id)
    {
        if (internal_id<0) return;

        // register variable into shared
        size_t var_id = _vars.registerVarialbe(var);
        // make a link between shared and system variable
        _map_sharedVar2sysVar.push_back(std::make_pair(var, std::make_pair(sys, internal_id)));
        // update a list of systems to be executed
        if (std::find(_listChildren.begin(), _listChildren.end(), sys) == _listChildren.end()) {
            _listChildren.push_back(sys);
        }
    }

    void add(const std::string &var)
    {
        // register variable into shared
        size_t var_id = _vars.registerVarialbe(var);
        _list_input_var.push_back(var_id);
    }

    int getVariableID(const std::string &var)
    {
        return _vars.getVariableID(var);
    }

    /// connect system input and shared variable
    void connectInput(const std::string &var, ICouplingProblem* sys, size_t internal_id)
    {
        //int varId = _vars.getVariableID(var);
        //if (varId<0) return;
        _map_sys2input[sys].push_back(make_pair(internal_id, var));
    }

    /// set input parameter 
    void set(size_t in_var, Variable* var)
    {
        _vars.setVariable(in_var, var);
    }

    /// get output variables
    Variable* get(size_t out) const
    {
        return _vars.getVariable(out);
    }

    size_t getNumberOfInputVarameters() const {return _list_input_var.size();};

    size_t getNumberOfOutputParameters() const {return _vars.size();};

    /// update shared variables from systems
    void updateSharedVariables()
    {
        // update own variables
        const size_t n_vars = _map_sharedVar2sysVar.size();
        for (size_t i=0; i<n_vars; ++i) {
            const std::string &var_name = _map_sharedVar2sysVar[i].first;
            PairSysVarId &x = _map_sharedVar2sysVar[i].second;
            ICouplingProblem *func = x.first;
            size_t f_var_id = x.second;
            // get v
            Variable* v = func->get(f_var_id);
            // share v
            _vars.setVariable(var_name, v);
        }

        // 
    }

    /// check consistency
    bool check()
    {
        bool flag = true;

        for (size_t i=0; i<_listChildren.size(); i++) {
            if (!_listChildren[i]->check()) 
                flag = false;
        }

        for (size_t i=0; i<_listChildren.size(); i++) {
            ICouplingProblem* problem = _listChildren[i];
            std::vector<PairInputVar> &vec = _map_sys2input[problem];
            for (size_t j=0; j<vec.size(); j++) {
                size_t internalId = vec[j].first;
                std::string &var_name = vec[j].second;
                if (!_vars.hasVariable(var_name)) {
                    std::cout << "*** Error: Variable '" << var_name << "' is not set." << std::endl;
                    flag = false;
                }
            }
        }

        return flag;
    }

    /// solve this system
    int solve()
    {
        return _algorithm->solve(_listChildren, _map_sharedVar2sysVar, &_vars);
    }

};
}
