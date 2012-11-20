/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceMid.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

//# Mid-plane elements

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceMidQuad(double* N4, const double* u)
{
    ShapeFunctionLine(N4, u);

    int i;
    //  N4[0] = N1[0];
    //  N4[1] = N1[1];
    N4[2] = N4[1];
    N4[3] = N4[0];
    for (i = 0; i < 4; i++)
        N4[i] *= 0.5;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceMidQuadHQ(double* N4, const double* u)
{
    ShapeFunctionLineHQ(N4, u);

    int i;
    //  N4[0] = N1[0];
    //  N4[1] = N1[1];
    N4[2] = N4[1];
    N4[3] = N4[0];
    for (i = 0; i < 4; i++)
        N4[i] *= 0.5;

    N4[0] = -N4[0];
    N4[1] = -N4[1];
    N4[4] = -N4[2];
    N4[3] = N4[0];
    N4[2] = N4[1];
    N4[5] = N4[2];
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceMidQuad(double* dN8, const double* u)
{
    //dN/dr
    GradShapeFunctionLine(dN8, u);
    //dN4[0] = dN4[0];
    //dN4[1] = dN4[1];
    dN8[2] = dN8[1];
    dN8[3] = dN8[0];
    for (int i = 0; i < 4; i++)
        dN8[i] *= 0.5;
    //dN/ds
    for (int i = 4; i < 8; i++)
        dN8[i] = 0.0;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceMidTri(double* N3, const double* u)
{
    int i;
    //*Notice: node order of an element
    //          2
    //     4 _ -*
    //   _ -*   |
    //  *---*---*
    //  0   3    1

    ShapeFunctionLine(N3, u);

    N3[2] = N3[1];
    for (i = 1; i < 3; i++)
        N3[i] *= 0.5;
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceMidTri(double* dN6, const double* u)
{
    // dN/dr
    GradShapeFunctionLine(dN6, u);
    dN6[2] = dN6[1];
    for (int i = 1; i < 3; i++)
        dN6[i] *= 0.5;
    // dN/ds
    for (int i = 3; i < 6; i++)
        dN6[i] = 0.0;
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2010 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceMidTriEnd(double* dN6, const double* u)
{
    //*Notice: node order of an element
    //  2
    //  *- _
    //  |    - _
    //  *-------*
    //  0       1
    // dN/dr
    double dN2[2] = {};
    GradShapeFunctionLine(dN2, u);
    dN6[0] = 0.5 * dN2[0];
    dN6[1] = dN2[1];
    dN6[2] = dN6[0];
    // dN/ds
    for (int i = 3; i < 6; i++)
        dN6[i] = 0.0;
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceMidPrism(double* N6, const double* u)
{
    double N3[3] = {};
    ShapeFunctionTri(N3, u);

    //*Notice: node order of an element
    // (bottom surface) (top surface)
    //      2             5
    //      *             *
    //     / \           / \
    //  8 *   * 7    11 *   * 10
    //   /     \       /     \
    //  *---*---*     *---*---*
    //  0   6    1    3   9   4
    for (int i = 0; i < 3; i++)
    {
        N6[i] = 0.5 * N3[i];
        N6[3 + i] = 0.5 * N3[i];
    }
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceMidPrism(double* dN12, const double* u)
{
    // dN/dr
    double dN6[3] = {};
    GradShapeFunctionTri(dN6, u);

    for (int i = 0; i < 3; i++)
    {
        dN12[i] = 0.5 * dN6[i];
        dN12[3 + i] = 0.5 * dN6[i];
        dN12[6 + i] = 0.5 * dN6[3 + i];
        dN12[9 + i] = 0.5 * dN6[3 + i];
    }
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceMidPyramid(double* N5, const double* u)
{
    double N3[3] = {};
    ShapeFunctionTri(N3, u);

    //*Notice: node order of an element
    // (top surface) (bottom surface)
    //      4             4
    //      *             *
    //     / \           / \
    //  7 *   * 6    10 *   * 9
    //   /     \       /     \
    //  *---*---*     *---*---*
    //  0   5    1    3   8    2
    N5[0] = 0.5 * N3[0];
    N5[1] = 0.5 * N3[1];
    N5[2] = 0.5 * N3[1];
    N5[3] = 0.5 * N3[0];
    N5[4] = N3[2];
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceMidPyramid(double* dN10, const double* u)
{
    // dN/dr
    double dN6[3] = {};
    GradShapeFunctionTri(dN6, u);

    for (int i = 0; i < 2; i++)
    {
        dN10[0 + 5 * i] = 0.5 * dN6[0 + 3 * i];
        dN10[1 + 5 * i] = 0.5 * dN6[1 + 3 * i];
        dN10[2 + 5 * i] = 0.5 * dN6[1 + 3 * i];
        dN10[3 + 5 * i] = 0.5 * dN6[0 + 3 * i];
        dN10[4 + 5 * i] = dN6[2 + 3 * i];
    }
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionInterfaceMidTet(double* N4, const double* u)
{
    double N3[3] = {};
    ShapeFunctionTri(N3, u);

    //*Notice: node order of an element
    // (top surface) (bottom surface)
    //      2             2
    //      *             *
    //     / \           / \
    //  6 *   * 5     8 *   * 5
    //   /     \       /     \
    //  *---*---*     *---*---*
    //  0   4    1    3   7    1
    N4[0] = 0.5 * N3[0];
    N4[1] = N3[1];
    N4[2] = N3[2];
    N4[3] = 0.5 * N3[0];
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionInterfaceMidTet(double* dN8, const double* u)
{
    // dN/dr
    double dN6[3] = {};
    GradShapeFunctionTri(dN6, u);

    for (int i = 0; i < 2; i++)
    {
        dN8[0 + 4 * i] = 0.5 * dN6[0 + 3 * i];
        dN8[1 + 4 * i] = dN6[1 + 3 * i];
        dN8[2 + 4 * i] = dN6[2 + 3 * i];
        dN8[3 + 4 * i] = 0.5 * dN6[0 + 3 * i];
    }
}



