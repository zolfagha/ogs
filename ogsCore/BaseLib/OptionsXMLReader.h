/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OptionsXMLReader.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

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
