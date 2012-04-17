
#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <cassert>

#include "Base/CodingTools.h"

namespace MathLib
{

class INamedIO
{
public:
	virtual ~INamedIO() {};

    virtual void setInputParameterName(size_t i, const std::string& key) = 0;

    virtual void setOutputParameterName(size_t i, const std::string& key) = 0;

    virtual bool isValid() const = 0;

    virtual bool hasInputParameter(const std::string& key) const = 0;

    virtual bool hasOutputParameter(const std::string& key) const = 0;

    /// get the name of parameter with the given name
    virtual const std::string& getInputParameterName(size_t i) const = 0;

    /// get the name of parameter with the given name
    virtual const std::string& getOutputParameterName(size_t i) const = 0;

    /// get the index of parameter with the given name
    virtual int getInputParameterID(const std::string&  key) const = 0;

    /// get the index of parameter with the given name
    virtual int getOutputParameterID(const std::string&  key) const = 0;

    /// get the number of output parameters
    virtual size_t getNumberOfInputParameterNames() const = 0;

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameterNames() const = 0;
};

template <typename T_BASE>
class NamedIO : public T_BASE
{
public:
	NamedIO() {};
	virtual ~NamedIO() {};

    virtual void setInputParameterName(size_t i, const std::string& key)
    {
    	if (i+1>_in_para_names.size()) _in_para_names.resize(i+1);
    	_in_para_names[i] = key;
    }

    virtual void setOutputParameterName(size_t i, const std::string& key)
    {
    	if (i+1>_out_para_names.size()) _out_para_names.resize(i+1);
    	_out_para_names[i] = key;
    }

    virtual bool isValid() const
    {
    	if (T_BASE::getNumberOfInputParameters()!=_in_para_names.size()) {
    		std::cout << "***Error: the number of input parameter names is not equal to that of values." << std::endl;
    		return false;
    	}
    	if (T_BASE::getNumberOfOutputParameters()!=_out_para_names.size()) {
    		std::cout << "***Error: the number of output parameter names is not equal to that of values." << std::endl;
    		return false;
    	}
    	return true;
    }


    virtual bool hasInputParameter(const std::string& key) const
    {
    	return _in_para_names.end() != std::find(_in_para_names.begin(), _in_para_names.end(), key);
    }

    virtual bool hasOutputParameter(const std::string& key) const
    {
    	return _out_para_names.end() != std::find(_out_para_names.begin(), _out_para_names.end(), key);
    }

    /// get the name of parameter with the given name
    virtual const std::string& getInputParameterName(size_t i) const { return _in_para_names[i];};

    /// get the name of parameter with the given name
    virtual const std::string& getOutputParameterName(size_t i) const { return _out_para_names[i];};

    /// get the index of parameter with the given name
    virtual int getInputParameterID(const std::string&  key) const
    {
    	std::vector<std::string>::const_iterator itr = std::find(_in_para_names.begin(), _in_para_names.end(), key);
    	if (itr==_in_para_names.end()) {
    		return -1;
    	} else {
    	   	return itr - _in_para_names.begin();
    	}
    };

    /// get the index of parameter with the given name
    virtual int getOutputParameterID(const std::string&  key) const
    {
    	std::vector<std::string>::const_iterator itr = std::find(_out_para_names.begin(), _out_para_names.end(), key);
    	if (itr==_out_para_names.end()) {
    		return -1;
    	} else {
    	   	return itr - _out_para_names.begin();
    	}
    }

    /// get the number of output parameters
    virtual size_t getNumberOfInputParameterNames() const {return _in_para_names.size();};

    /// get the number of output parameters
    virtual size_t getNumberOfOutputParameterNames() const {return _out_para_names.size();};

private:
    std::vector<std::string> _in_para_names;
    std::vector<std::string> _out_para_names;
};
} //end
