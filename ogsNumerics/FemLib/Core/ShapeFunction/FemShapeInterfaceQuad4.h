/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceQuad4.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"
#include "FemShapeLine2.h"

namespace FemLib
{

/**
 * Goodman joint elements: [[u]] = u_top - u_bottom
 *
 *  3                2
 *  *----------------*
 *  |                |
 *  *----------------*
 *  0                1
 */
class FemShapeInterfaceQuad4 : public TemplateShapeFunction<2, 4>
{
public:
    void computeShapeFunction(const double* pt, double *N4)
    {
        double N2[2] = {};
        _fe_line2.computeShapeFunction(pt, N2);

        //*Notice: node order of an element
        //  3                2
        //  *----------------*
        //  |                |
        //  *----------------*
        //  0                1
        N4[0] = -N2[0];
        N4[1] = -N2[1];
        N4[3] = N2[0];
        N4[2] = N2[1];
    };

    void computeGradShapeFunction(const double*, double *dN)
    {
        //TODO
        dN[0] = -0.5;
        dN[1] = 0.5;
    };

private:
    FemShapeLine2 _fe_line2;
};

}

