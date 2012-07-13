
#pragma once

#include <string>
#include <map>

#include "BaseLib/CodingTools.h"
//#include "BaseLib/Options.h"
#include "ProcessInfo.h"
#include "ProcessFactoryBase.h"
#include "ProcessFactoryImpl.h"
#include "Process.h"


namespace ProcessLib
{


/**
 * Singleton class
 */
class AbstractProcessBuilder
{

protected:
//#define PROCRSS_REGISTER
    AbstractProcessBuilder() {};
public:
	virtual ~AbstractProcessBuilder() {};

	ProcessInfo* registerProcess(const std::string &pcs_name, ProcessFactoryBase* pcs_buid);

	bool hasRegisterd(const std::string &pcs_name) const;

	Process* create(const std::string &pcs_name) const;

	void output() const;


private:
	std::map<std::string, ProcessFactoryBase*> _map_pcs_name2new;
};

} //end

#ifdef PROCRSS_REGISTER

#define OGS_PROCESS(pcs_name, pcs_classname)\
ProcessLib::ProcessInfo* const pcs_classname::_pcs_info =\
		ProcessLib::ProcessBuilder::getInstance()->registerProcess(\
        #pcs_name,\
        new ProcessLib::ProcessFactoryImpl<pcs_classname >);

#define OGS_DEF_PROCESS(pcs_name, pcs_classname)\
int class_##pcs_name = 0;

#define OGS_LINK_PROCESS(pcs_name, pcs_classname)\
extern int class_##pcs_name; \
int link_##pcs_name = class_##pcs_name;

#else

#define OGS_PROCESS(pcs_name, pcs_classname)
#define OGS_DEF_PROCESS(pcs_name, pcs_classname)
#define OGS_LINK_PROCESS(pcs_name, pcs_classname)
#define OGS_ADD_PROCESS(pcs_name, pcs_classname)\
		registerProcess(#pcs_name, new ProcessLib::ProcessFactoryImpl<pcs_classname>);

#endif

