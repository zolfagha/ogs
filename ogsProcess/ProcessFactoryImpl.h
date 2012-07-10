
#pragma once

#include <string>

#include "ProcessFactoryBase.h"

namespace ProcessLib
{

// This class provides implementation of TeastFactoryBase interface.
// It is used in TEST and TEST_F macros.
template <class T>
class ProcessFactoryImpl : public ProcessFactoryBase
{
public:
  virtual Process* createProcess() { return new T; }
};

}
