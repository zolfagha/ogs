
#pragma once

#include "BaseLib/CodingTools.h"

namespace ProcessLib
{

// forward declaration
class Process;

/**
 * \brief Defines the abstract factory interface that creates instances
 * of a Test object.
 */
class ProcessFactoryBase
{
public:
  virtual ~ProcessFactoryBase() {}

  /// Creates a test instance to run. The instance is both created and destroyed
  /// within TestInfoImpl::Run()
  virtual Process* createProcess() = 0;

protected:
  ProcessFactoryBase() {}

private:
  DISALLOW_COPY_AND_ASSIGN(ProcessFactoryBase);
};

}
