/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessInfo.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */


#pragma once

#include <string>

namespace ProcessLib
{

/**
 *
 */
class ProcessInfo
{
public:
    ~ProcessInfo() {};

    const char* getName() const { return _name.c_str(); }

private:
    std::string _name;
};

}
