/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFemNeumannBC.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>


namespace SolutionLib
{

/**
 *
 */
class IFemNeumannBC
{
public:
    ///
    virtual ~IFemNeumannBC() {};

    virtual void setOrder(size_t order) = 0;

    /// setup B.C.
    virtual void setup() = 0;

    ///
    virtual std::vector<size_t>& getListOfBCNodes() = 0;

    ///
    virtual std::vector<double>& getListOfBCValues() = 0;
};

}
