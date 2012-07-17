
#pragma once

#include "Process.h"
#include "ProcessFactoryBase.h"

namespace ProcessLib
{

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
