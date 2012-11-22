/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceTri3.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

#include "TemplateShapeFunction.h"

namespace FemLib
{


class FemShapeInterfaceTri3 : public TemplateShapeFunction<2, 3>
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = 1.0 - pt[0];
        N[1] = 1.0 + pt[0];
        for (int i = 0; i < 2; i++)
            N[i] *= 0.5;
    };

    void computeGradShapeFunction(const double*, double *dN)
    {
        dN[0] = -0.5;
        dN[1] = 0.5;
    };
};

}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   05/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionGoodmanTri(double* N3, const double* u)
{
    double N2[2];
    ShapeFunctionLine(N2, u);

    //*Notice: node order of an element
    //          2
    //       _ -*
    //   _ -    |
    //  *-------*
    //  0       1

    N3[0] = 0.0;
    N3[1] = -N2[1];
    N3[2] = N2[1];
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   05/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionGoodmanTriEnd(double* N3, const double* u)
{
    double N2[2] = {};
    ShapeFunctionLine(N2, u);

    //*Notice: node order of an element
    //  2
    //  *- _
    //  |    - _
    //  *-------*
    //  0       1

    N3[0] = -N2[1];
    N3[1] = 0.0;
    N3[2] = N2[1];
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   05/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionGoodmanTriHQ(double* N5, const double* u)
{
    double N3[3];
    ShapeFunctionLineHQ(N3, u);

    //*Notice: node order of an elemet
    //          2
    //     4 _ -*
    //   _ -*   |
    //  *---*---*
    //  0   3    1

    N5[0] = 0.0;
    N5[1] = -N3[1];
    N5[3] = -N3[2];
    N5[2] = N3[1];
    N5[4] = N3[2];
}

