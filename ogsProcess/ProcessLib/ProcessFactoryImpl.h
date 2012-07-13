
#pragma once

#include "ProcessFactoryBase.h"

namespace ProcessLib
{

// forward declaration
class Process;

/**
 * \brief This class provides implementation of TeastFactoryBase interface.
 * It is used in TEST and TEST_F macros.
 */
template <class T>
class ProcessFactoryImpl : public ProcessFactoryBase
{
public:
  virtual Process* createProcess() { return new T; }
};

}
