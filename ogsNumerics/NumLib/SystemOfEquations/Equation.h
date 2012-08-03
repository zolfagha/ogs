/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Equation.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <valarray>

#include "BaseLib/CodingTools.h"
#include "NumLib/DataType.h"
#include "Variable.h"


namespace NumLib
{

struct TimeODEComponent
{
    LocalMatrix* mass;
    LocalMatrix* stiffness;
    LocalVector* source;
};

struct TimeODE
{
    std::vector<Variable> variables;
    std::vector<TimeODEComponent> componets;
    Variable primary_variable;
    size_t n_dof;

    void addVariable(Variable &var, TimeODEComponent &comp, bool is_primary=false)
    {
        variables.push_back(var);
        componets.push_back(comp);
        if (is_primary) {
            primary_variable = var;
            n_dof = var.n_dof;
        }
    }

};

struct LinearComponent
{
    LocalMatrix* stiffness;
    LocalVector* source;

    LinearComponent(LocalMatrix* m, LocalVector* v = 0) : stiffness(m), source(v) {};
};

struct LinearEquation
{
    std::vector<Variable*> variables;
    std::vector<LinearComponent*> componets;
    Variable* primary_variable;
    size_t n_dof;

    void addVariable(Variable &var, LinearComponent &comp, bool is_primary=false)
    {
        variables.push_back(&var);
        componets.push_back(&comp);
        if (is_primary) {
            primary_variable = &var;
            n_dof = var.n_dof;
        }
    }

    virtual ~LinearEquation()
    {
        BaseLib::releaseObjectsInStdVector(componets);
    }
};

} //end

