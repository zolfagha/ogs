
#pragma once

#include <string>
#include <map>

#include "BaseLib/CodingTools.h"

namespace NumLib
{
class TemplateTransientMonolithicSystem;
}


namespace ProcessLib
{

typedef NumLib::TemplateTransientMonolithicSystem Process;

//================================================================================
// Copy and modify from gtest-internal.h
//================================================================================
class ProcessInfo
{
public:
	~ProcessInfo() {};
    const char* name() const { return _name.c_str(); }

private:
    std::string _name;
};

// Defines the abstract factory interface that creates instances
// of a Test object.
class ProcessFactoryBase
{
public:
  virtual ~ProcessFactoryBase() {}

  // Creates a test instance to run. The instance is both created and destroyed
  // within TestInfoImpl::Run()
  virtual Process* createProcess() = 0;

protected:
  ProcessFactoryBase() {}

private:
  DISALLOW_COPY_AND_ASSIGN(ProcessFactoryBase);
};

// This class provides implementation of TeastFactoryBase interface.
// It is used in TEST and TEST_F macros.
template <class ProcessClass>
class ProcessFactoryImpl : public ProcessFactoryBase
{
public:
  virtual Process* createProcess() { return new ProcessClass; }
};

#define OGS_PROCESS(pcs_name, pcs_classname)\
ProcessLib::ProcessInfo* const pcs_classname::_pcs_info =\
   ProcessLib::ProcessBuilder::getInstance()->registerProcess(\
        #pcs_name,\
        new ProcessLib::ProcessFactoryImpl<pcs_classname >);\

//================================================================================

/**
 * Singleton class
 */
class ProcessBuilder
{
public:
	/// get a global instance
	static ProcessBuilder* getInstance();

	virtual ~ProcessBuilder() {};

	ProcessInfo* registerProcess(const std::string &pcs_name, ProcessFactoryBase* pcs_buid);

	bool hasRegisterd(const std::string &pcs_name) const;

	Process* create(const std::string &pcs_name) const;

	void output() const;

private:
    ProcessBuilder() {};

private:
    static ProcessBuilder* _obj;
	std::map<std::string, ProcessFactoryBase*> _map_pcs_name2new;
};


} //end
