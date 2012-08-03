/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ForwardEuler.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#pragma once

namespace MathLib
{

class ForwardEuler
{
public:
    void solve(double (*f)(double, double), double *t, int n, double *y) {
        for (int i=0; i<n-1; i++) {
            double dt = t[i+1] - t[i];
            y[i+1] = solve(f, y[i], t[i], dt);
        }
    };

    double solve(double (*f)(double, double), double y0, double t0, double dt) {
        double y = y0 + dt*f(t0, y0);
        return y;
    };
};
}
