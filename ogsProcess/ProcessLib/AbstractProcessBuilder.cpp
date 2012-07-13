
#include "AbstractProcessBuilder.h"

#include <iostream>
#include "logog/include/logog.hpp"

namespace ProcessLib
{

ProcessInfo* AbstractProcessBuilder::registerProcess(const std::string &pcs_name, ProcessFactoryBase* pcs_buid)
{
	this->_map_pcs_name2new[pcs_name] = pcs_buid;
	ProcessInfo* pcs_info = new ProcessInfo();

	return pcs_info;
}

bool AbstractProcessBuilder::hasRegisterd(const std::string &pcs_name) const
{
	return _map_pcs_name2new.count(pcs_name)>0;
}

Process* AbstractProcessBuilder::create(const std::string &pcs_name) const
{
	if (!hasRegisterd(pcs_name)) return 0;

	std::map<std::string, ProcessFactoryBase*>::const_iterator itr = this->_map_pcs_name2new.find(pcs_name);
	return itr->second->createProcess();
}

void AbstractProcessBuilder::output() const
{
	std::map<std::string, ProcessFactoryBase*>::const_iterator itr;
    LOGOG_COUT << _LG("List of available modules") << std::endl;
	for (itr=this->_map_pcs_name2new.begin(); itr!=_map_pcs_name2new.end(); ++itr) {
        LOGOG_COUT << "* " << itr->first << std::endl;
	}
    LOGOG_COUT << std::endl << std::endl;
}

} //end
