
#pragma once

#include <string>
#include <map>
#include <vector>

#include "StringTools.h"
#include "BaseLib/Options.h"

namespace BaseLib
{

bool addXMLtoOptions(const std::string &xml_file, BaseLib::Options &properties);

}
