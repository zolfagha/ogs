/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessInfo.h
 *
 * Created on 2010-09-02 by Thomas Fischer
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
