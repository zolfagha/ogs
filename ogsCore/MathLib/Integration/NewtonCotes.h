
#pragma once

namespace MathLib
{

class NewtonCotes
{
    double middpoint(double (*fun)(double), double a, double b, size_t n) {
        double h = (b-a) / (double)n;
        double val = .0;
        double x = a + 0.5*h;
        for (int i=0; i<n; i++) {
            val += fun(x);
            x += h;
        }
        val *= h;
        return val;
    };

    double trapezoid(double (*fun)(double), double a, double b, size_t n) {
        double h = (b-a) / (double)n;
        double val = 0.5*(fun(a) + fun(b));
        double x = a;
        for (int i=0; i<n-1; i++) {
            x += h;
            val += fun(x);
        }
        val *= h;
        return val;
    };

    double simpson(double (*fun)(double), double a, double b, size_t n) {
        double h = (b-a) / (double)n;
        double val = fun(a) + 4*fun(b-0.5*h) + fun(b);
        double x = a;
        for (int i=0; i<n-1; i++) {
            val += 4*fun(x+0.5*h) + 2*fun(x+h);
            x += h;
        }
        val = val * h / 6.0;
        return val;
    };
};

}
