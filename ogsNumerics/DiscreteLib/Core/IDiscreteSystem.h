/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace DiscreteLib
{

/**
 * \brief Interface for all kinds of discrete systems
 *
 *  Discrete system class contains the followings
 *  - discrete points (space, time) 
 *  - discrete data (i.e. vector)
 *  - linear equations which are used to calculate discrete data
 */
class IDiscreteSystem
{
public:
    virtual ~IDiscreteSystem() {};
    
    //IDiscreteLinearEquation* createLinearEquation();
    //IDiscreteVector* createVector();
};

} //end
