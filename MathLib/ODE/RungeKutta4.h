
#pragma once

namespace MathLib
{

// dy/dt = f(t, y)
class RungeKutta4
{
public:
    static void solve(double (*f)(double, double), double *t, int n, double *y) {
        for (int i=0; i<n-1; i++) {
            double dt = t[i+1] - t[i];
            y[i+1] = solve(f, y[i], t[i], dt);
        }
    };

    static double solve(double (*f)(double, double), double y0, double t0, double dt) {

        double k1 = dt*f(t0,y0);
        double k2 = dt*f(t0+dt/2.0, y0+k1/2.0);
        double k3 = dt*f(t0+dt/2.0, y0+k2/2.0);
        double k4 = dt*f(t0+dt, y0+k3);
        double y = y0 + (k1+2.0*k2+2.0*k3+k4)/6.0;

        return y;
    }

};

}
