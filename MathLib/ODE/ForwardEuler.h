
#pragma once

namespace MathLib
{

class ForwardEuler
{
public:
    static void solve(double (*f)(double, double), double *t, int n, double *y) {
        for (int i=0; i<n-1; i++) {
            double dt = t[i+1] - t[i];
            y[i+1] = solve(f, y[i], t[i], dt);
        }
    };

    static double solve(double (*f)(double, double), double y0, double t0, double dt) {
        double y = y0 + dt*f(t0, y0);
        return y;
    };
};
}
