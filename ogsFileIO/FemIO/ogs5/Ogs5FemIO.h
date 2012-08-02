
#pragma once

#include <string>
#include <vector>

#include "Ogs5FemData.h"

namespace ogs5
{

namespace Ogs5FemIO
{
bool read(const std::string &proj_path, Ogs5FemData &ogs5data);
};

}
