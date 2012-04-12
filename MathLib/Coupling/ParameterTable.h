
#pragma once

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>

#include "Base/CodingTools.h"
#include "Parameter.h"

namespace MathLib
{

/**
 * \brief Container of named variables
 */
class ParameterSet
{
public:
	///
	ParameterSet() {};

    ///
    virtual ~ParameterSet()
    {
    	clear();
    }

    /// reset data
    void clear()
    {
        _list_var_data.clear();
        _list_attributes.clear();
        Base::releaseObjectsInStdVector(_list_own_data);
    }

    ///// make a copy of the given object
    //void assign(const ParameterSet &src)
    //{
    //	this->clear();
    //    _list_var_names.assign(src._list_var_names.begin(), src._list_var_names.end());
    //    _list_var_data.resize(src._list_var_data.size());
    //    for (size_t i=0; i<src._list_var_data.size(); i++) {
    //        if (src._list_var_data[i]!=0)  {
    //            _list_var_data[i] = src._list_var_data[i]->clone();
    //        }
    //    }
    //}

    void move(ParameterSet &dest)
    {
        dest.clear();
        dest._list_attributes.assign(this->_list_attributes.begin(), this->_list_attributes.end());
        //if (copy_val) {
        //    dest._list_var_data.resize(this->_list_var_data.size());
        //    for (size_t i=0; i<this->_list_var_data.size(); i++) {
        //        if (this->_list_var_data[i]!=0)  {
        //            dest._list_var_data[i] = this->_list_var_data[i]->clone();
        //        }
        //    }
        //} else {
            dest._list_var_data.assign(this->_list_var_data.begin(), this->_list_var_data.end());
        //}
        // move ownership of memory
        dest._list_own_data.assign(this->_list_own_data.begin(), this->_list_own_data.end());
        _list_own_data.clear();
    }

    /// return the size of the stores variables
    size_t size() const 
    {
        return _list_var_data.size();
    }

    /// return if the container has a variable with the given key
    bool contain(const std::string &var_name) const
    {
    	return (this->find(var_name)>=0);
    }

    /// return the name of the variable with the given index
    const std::string& getName(size_t var_id) const 
    {
        return _list_attributes[var_id].name;
    };

    /// return the index of the variable with the given name
    int find(const std::string& name) const 
    {
        for (std::vector<ParameterAttribute>::const_iterator itr = _list_attributes.begin(); itr!=_list_attributes.end(); ++itr) {
            if (itr->name.compare(name)==0) {
                return (itr - _list_attributes.begin());
            }
        }

        return -1;
    }

    /// return the variable with the given index
    const Parameter* get(size_t var_id) const
    {
        return _list_var_data[var_id];
    }

    /// return the variable with the given index
    template <class T>
    const T* get(size_t var_id) const
    {
        return static_cast<const T*>(_list_var_data[var_id]);
    }

    /// register the variable
    size_t add(const std::string &var_name)
    {
        size_t new_id = _list_attributes.size();
        _list_attributes.resize(_list_attributes.size()+1);
        _list_attributes[new_id].name = var_name;
        _list_var_data.push_back(0);
        return new_id;
    }

    size_t addInput(const std::string &var_name)
    {
        size_t new_id = _list_attributes.size();
        _list_attributes.resize(_list_attributes.size()+1);
        _list_attributes[new_id].name = var_name;
        _list_attributes[new_id].is_fxied = true;
        _list_var_data.push_back(0);
        return new_id;
    }

    /// set variable
    void set(size_t var_id,  const Parameter& v)
    {
        //assert(&v!=0);
        if (&v==0) return;
        if (isFixed(var_id)) {
            _list_var_data[var_id] = &v;
        } else {
            _list_var_data[var_id] = v.clone();
            _list_own_data.push_back(_list_var_data[var_id]);
        }
        _list_attributes[var_id].is_updated = true;
    }

    bool isFixed(size_t var_id) const
    {
        return _list_attributes[var_id].is_fxied;
    }

    void setFixed(size_t var_id, bool is_fxied)
    {
        _list_attributes[var_id].is_fxied = is_fxied;
    } 

    bool isUpdated(size_t var_id) const
    {
        return _list_attributes[var_id].is_updated;
    }

    void finishUpdate()
    {
        for (size_t i=0; i<_list_attributes.size(); i++)
            _list_attributes[i].is_updated = false;
    } 

private:
    struct ParameterAttribute
    {
        std::string name;
        bool is_fxied;
        bool is_updated;
        ParameterAttribute() : is_fxied(false), is_updated(false) {};
    };
    std::vector<const Parameter*> _list_var_data;
    std::vector<ParameterAttribute> _list_attributes;
    std::vector<const Parameter*> _list_own_data;

    DISALLOW_COPY_AND_ASSIGN(ParameterSet);
};

}
