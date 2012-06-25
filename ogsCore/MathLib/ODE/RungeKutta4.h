
#pragma once

namespace MathLib
{

// dy/dt = f(t, y)
template <class F_DXDT>
class RungeKutta4
{
public:
	template <class T_TIM_VECTOR, class T_Y_VECTOR>
    void solve(F_DXDT &f, const T_TIM_VECTOR &t, int n, T_Y_VECTOR &y)
    {
        for (int i=0; i<n-1; i++) {
            double dt = t[i+1] - t[i];
            solve(f, t[i], dt, y[i], y[i+1]);
        }
    };

	template <class T_VALUE>
    void solve(F_DXDT &f, double t0, double dt, const T_VALUE &y0, T_VALUE &y)
    {
		T_VALUE k1 = dt*f(t0,y0);
		T_VALUE k2 = dt*f(t0+dt/2.0, y0+k1/2.0);
		T_VALUE k3 = dt*f(t0+dt/2.0, y0+k2/2.0);
		T_VALUE k4 = dt*f(t0+dt, y0+k3);
        y = y0 + (k1+2.0*k2+2.0*k3+k4)/6.0;
    }
};

} //end
