
#pragma once

#include <string>
#include <map>

namespace ProcessLib
{

class ProcessBuilder
{
private:
	static ProcessBuilder* _obj;
	ProcessBuilder() {};
public:
	static ProcessBuilder* getObject()
	{
		if (_obj==0) _obj = new ProcessBuilder();
		return _obj;
	}

	virtual ~ProcessBuilder() {};

	void registerProcess(const std::string &pcs_name, void* F_CREATE_PCS)
	{
		this->_map_pcs_name2new[pcs_name] = F_CREATE_PCS;
	}

	bool hasRegisterd(const std::string &pcs_name) const
	{
		return _map_pcs_name2new.count(pcs_name)>0;
	}

	void* create(const std::string &pcs_name) const
	{
		if (!hasRegisterd(pcs_name)) return 0;

		std::map<std::string, void*>::const_iterator itr = this->_map_pcs_name2new.find(pcs_name);
		return itr->second;
	}

private:
	std::map<std::string, void*> _map_pcs_name2new;
};

} //end
