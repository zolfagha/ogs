/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Process.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "logog.hpp"

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"

namespace ProcessLib
{

/**
 * \brief Interface definition of process class
 *
 */
//typedef NumLib::AbstractTransientMonolithicSystem Process;

/**
 * \brief Interface definition of process class
 *
 */
class Process : public NumLib::AbstractTransientMonolithicSystem
{
public:
    Process(const std::string &pcs_name, size_t n_in_parameters, size_t n_out_parameters)
        : _pcs_name(pcs_name)
    {
        resizeInputParameter(n_in_parameters);
        resizeOutputParameter(n_out_parameters);
    }

    virtual ~Process() {};

    const std::string& getProcessName() const {return _pcs_name;};

private:
    std::string _pcs_name;
};


}

