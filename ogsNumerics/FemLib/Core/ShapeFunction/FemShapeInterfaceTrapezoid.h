/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceTrapezoid.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

//# Trapezoidal elements

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceTrapezoidQuad(double* N8, const double* u)
{
    //*Notice: node order of an element
    //  3                2
    //  *----------------*
    //  |                |
    //  *----------------*
    //  0                1

    double N2[3] = {};
    ShapeFunctionLine(N2, u);

    N8[0] = 0.0;
    N8[1] = 0.0;
    N8[2] = N2[1];
    N8[3] = N2[0];
    N8[4] = N2[0];
    N8[5] = N2[1];
    N8[6] = 0.0;
    N8[7] = 0.0;
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceTrapezoidQuad(double* dN8, const double* u)
{
    //*Notice: node order of an element
    //  3                2
    //  *----------------*
    //  |                |
    //  *----------------*
    //  0                1

    double dN2[2] = {};
    GradShapeFunctionLine(dN2, u);

    dN8[0] = 0.0;
    dN8[1] = 0.0;
    dN8[2] = dN2[1];
    dN8[3] = dN2[0];
    dN8[4] = dN2[0];
    dN8[5] = dN2[1];
    dN8[6] = 0.0;
    dN8[7] = 0.0;
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceTrapezoidTri(double* N6, const double* u)
{
    //*Notice: node order of an element
    //          2
    //       _ -*
    //   _ -    |
    //  *-------*
    //  0        1

    double N2[3] = {};
    ShapeFunctionLine(N2, u);

    N6[0] = N2[0];
    N6[1] = 0.0;
    N6[2] = N2[1];
    N6[3] = N2[0];
    N6[4] = N2[1];
    N6[5] = 0.0;
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceTrapezoidTri(double* dN6, const double* u)
{
    //*Notice: node order of an element
    //          2
    //       _ -*
    //   _ -    |
    //  *-------*
    //  0        1

    double dN2[2] = {};
    GradShapeFunctionLine(dN2, u);

    dN6[0] = dN2[0];
    dN6[1] = 0.0;
    dN6[2] = dN2[1];
    dN6[3] = dN2[0];
    dN6[4] = dN2[1];
    dN6[5] = 0.0;
}


