/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfacePrism12.h
 *
 * Created on 2012-11-20 by Norihiro Watanabe
 */

#pragma once

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   05/2010 NW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionGoodmanPrismHQ(double* N12, const double* u)
{
    double N6[6] = {};
    ShapeFunctionTriHQ(N6, u);

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
        N12[i] = -N6[i];
        N12[6 + i] = -N6[3 + i];
        N12[3 + i] = N6[i];
        N12[9 + i] = N6[3 + i];
    }
}
