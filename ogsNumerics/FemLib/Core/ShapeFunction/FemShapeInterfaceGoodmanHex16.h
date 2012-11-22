/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemShapeInterfaceGoodmanHex16.h
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
void ShapeFunctionGoodmanHexHQ(double* N16, const double* u)
{
    double N8[8] = {};
    ShapeFunctionQuadHQ8(N8, u);

    // Node order of an element
    //      (top surface)         (bottom surface)
    //  3        10       2      7        14       6
    //  *--------*--------*      *--------*--------*
    //  |                 |      |                 |
    //  * 11              * 9    * 15              * 13
    //  |                 |      |                 |
    //  *--------*--------*      *--------*--------*
    //  0        8        1      4        12       5
    for (int i = 0; i < 4; i++)
    {
        N16[i] = N8[i];
        N16[8 + i] = N8[4 + i];
        N16[4 + i] = -N8[i];
        N16[12 + i] = -N8[4 + i];
    }
}

