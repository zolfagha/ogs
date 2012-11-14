
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "TemplateShapeFunction.h"

namespace FemLib
{


class FemShapeLine2 : public TemplateShapeFunction<1, 2> 
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

class FemShapeLine3 : public TemplateShapeFunction<1, 3> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = 0.5 * pt[0] * (pt[0] - 1.0);
        N[1] = 0.5 * pt[0] * (pt[0] + 1.0);
        N[2] = 1.0 - pt[0] * pt[0];
    };

    void computeGradShapeFunction(const double* pt, double *dN)
    {
        dN[0] = pt[0] - 0.5;
        dN[1] = pt[0] + 0.5;
        dN[2] = -2.0 * pt[0];
    };
};

class FemShapeQuad4 : public TemplateShapeFunction<2, 4> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = (1.0 + pt[0]) * (1.0 + pt[1]);
        N[1] = (1.0 - pt[0]) * (1.0 + pt[1]);
        N[2] = (1.0 - pt[0]) * (1.0 - pt[1]);
        N[3] = (1.0 + pt[0]) * (1.0 - pt[1]);
        for (int i = 0; i < 4; i++)
            N[i] *= 0.25;
    };

    void computeGradShapeFunction(const double* u, double *dN4)
    {
        dN4[0] = +(1.0 + u[1]);
        dN4[1] = -(1.0 + u[1]);
        dN4[2] = -(1.0 - u[1]);
        dN4[3] = +(1.0 - u[1]);
        dN4[4] = +(1.0 + u[0]);
        dN4[5] = +(1.0 - u[0]);
        dN4[6] = -(1.0 - u[0]);
        dN4[7] = -(1.0 + u[0]);
        for (int i = 0; i < 8; i++)
            dN4[i] *= 0.25;
    };
};

class FemShapeQuad8 : public TemplateShapeFunction<2, 8> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[0] = -0.25 * (1.0 - pt[0]) * (1.0 - pt[1]) * (( 1.0 + pt[0] + pt[1]));
        N[1] =  0.25 * (1.0 + pt[0]) * (1.0 - pt[1]) * ((-1.0 + pt[0] - pt[1]));
        N[2] =  0.25 * (1.0 + pt[0]) * (1.0 + pt[1]) * ((-1.0 + pt[0] + pt[1]));
        N[3] = -0.25 * (1.0 - pt[0]) * (1.0 + pt[1]) * (( 1.0 + pt[0] - pt[1]));
        //
        N[4] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 - pt[1]);
        N[5] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 + pt[0]);
        N[6] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 + pt[1]);
        N[7] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 - pt[0]);
    };

    void computeGradShapeFunction(const double* pt, double *dN8)
    {
        double r = pt[0];
        double s = pt[1];

        //dN/dr
        dN8[0] = (1 - s) * (2 * r + s) * 0.25;
        dN8[1] = (1 - s) * (2 * r - s) * 0.25;
        dN8[2] = (1 + s) * (2 * r + s) * 0.25;
        dN8[3] = (1 + s) * (2 * r - s) * 0.25;
        dN8[4] = -r * (1 - s);
        dN8[5] = (1 - s * s) * 0.5;
        dN8[6] = -r * (1 + s);
        dN8[7] = -(1 - s * s) * 0.5;

        //dN/ds
        dN8[8] = (1 - r) * (r + 2 * s) * 0.25;
        dN8[9] = -(1 + r) * (r - 2 * s) * 0.25;
        dN8[10] = (1 + r) * (r + 2 * s) * 0.25;
        dN8[11] = -(1 - r) * (r - 2 * s) * 0.25;
        dN8[12] = -(1 - r * r) * 0.5;
        dN8[13] = -(1 + r) * s;
        dN8[14] = (1 - r * r) * 0.5;
        dN8[15] = -(1 - r) * s;
    };
};

class FemShapeQuad9 : public TemplateShapeFunction<2, 9> 
{
public:
    void computeShapeFunction(const double* pt, double *N)
    {
        N[8] = (1.0 - pt[0] * pt[0]) * ( 1.0 - pt[1] * pt[1]);
        N[7] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 + pt[0]) - 0.5 * N[8];
        N[6] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 - pt[1]) - 0.5 * N[8];
        N[5] = 0.5 * (1.0 - pt[1] * pt[1]) * (1.0 - pt[0]) - 0.5 * N[8];
        N[4] = 0.5 * (1.0 - pt[0] * pt[0]) * (1.0 + pt[1]) - 0.5 * N[8];
        N[3] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]) - 0.5 * N[6] - 0.5 * N[7] - 0.25 * N[8];
        N[2] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]) - 0.5 * N[5] - 0.5 * N[6] - 0.25 * N[8];
        N[1] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]) - 0.5 * N[4] - 0.5 * N[5] - 0.25 * N[8];
        N[0] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]) - 0.5 * N[4] - 0.5 * N[7] - 0.25 * N[8];
    };

    void computeGradShapeFunction(const double* u, double *dN9)
    {
        dN9[8] = -2.0 * u[0] * (1.0 - u[1] * u[1]);
        dN9[7] = +0.5 * (1.0 - u[1] * u[1]) - 0.5 * dN9[8];
        dN9[6] = -1.0 * u[0] * (1.0 - u[1]) - 0.5 * dN9[8];
        dN9[5] = -0.5 * (1.0 - u[1] * u[1]) - 0.5 * dN9[8];
        dN9[4] = -1.0 * u[0] * (1.0 + u[1]) - 0.5 * dN9[8];
        dN9[3] = +0.25 * (1 - u[1]) - 0.5 * dN9[6] - 0.5 * dN9[7] - 0.25 * dN9[8];
        dN9[2] = -0.25 * (1 - u[1]) - 0.5 * dN9[5] - 0.5 * dN9[6] - 0.25 * dN9[8];
        dN9[1] = -0.25 * (1 + u[1]) - 0.5 * dN9[4] - 0.5 * dN9[5] - 0.25 * dN9[8];
        dN9[0] = +0.25 * (1 + u[1]) - 0.5 * dN9[4] - 0.5 * dN9[7] - 0.25 * dN9[8];

        dN9[17] = -2.0 * u[1] * (1.0 - u[0] * u[0]);
        dN9[16] = -1.0 * u[1] * (1.0 + u[0]) - 0.5 * dN9[17];
        dN9[15] = -0.5 * (1.0 - u[0] * u[0]) - 0.5 * dN9[17];
        dN9[14] = -1.0 * u[1] * (1.0 - u[0]) - 0.5 * dN9[17];
        dN9[13] = +0.5 * (1 - u[0] * u[0]) - 0.5 * dN9[17];
        dN9[12] = -0.25 * (1 + u[0]) - 0.5 * dN9[15] - 0.5 * dN9[16] - 0.25 * dN9[17];
        dN9[11] = -0.25 * (1 - u[0]) - 0.5 * dN9[14] - 0.5 * dN9[15] - 0.25 * dN9[17];
        dN9[10] = +0.25 * (1 - u[0]) - 0.5 * dN9[13] - 0.5 * dN9[14] - 0.25 * dN9[17];
        dN9[9] = +0.25 * (1 + u[0]) - 0.5 * dN9[13] - 0.5 * dN9[16] - 0.25 * dN9[17];
    };
};

class FemShapeHex8 : public TemplateShapeFunction<3, 8> 
{
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
};

class FemShapeHex20 : public TemplateShapeFunction<3, 20> 
{
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
};


class FemShapeTriangle3 : public TemplateShapeFunction<2, 3> 
{
    void computeShapeFunction(const double* u, double *N)
    {
        N[0] = 1. - u[0] - u[1];
        N[1] = u[0];
        N[2] = u[1];
    };
    void computeGradShapeFunction(const double*, double *dN3)
    {
        //   d()/dL_1
        dN3[0] = -1.0;
        dN3[1] =  1.0;
        dN3[2] =  0.0;
        //   d()/dL_2
        dN3[3] = -1.0;
        dN3[4] = 0.0;
        dN3[5] = 1.0;
    }
};

class FemShapeTriangle6 : public TemplateShapeFunction<2, 6> 
{
    void computeShapeFunction(const double* u, double *N6)
    {
        N6[0] = 2. * (1. - u[0] - u[1]) * (0.5 - u[0] - u[1]);
        N6[1] = u[0] * (2. * u[0] - 1.);
        N6[2] = u[1] * (2. * u[1] - 1.);
        N6[3] = 4. * u[0] * (1. - u[0] - u[1]);
        N6[4] = 4. * u[0] * u[1];
        N6[5] = 4. * u[1] * (1. - u[0] - u[1]);
    };
    void computeGradShapeFunction(const double* u, double *dN6)
    {
        dN6[0] = 4. * (u[0] + u[1]) - 3.;     // dN1/dL1
        dN6[6] = 4. * (u[0] + u[1]) - 3.;     // dN1/dL2

        dN6[1] = 4. * u[0] - 1.;              // dN2/dL1
        dN6[7] = 0.;                          // dN2/dL2

        dN6[2] = 0.;                          // dN3/dL1
        dN6[8] = 4. * u[1] - 1.;              // dN3/dL2

        dN6[3] =  4. * (1 - 2. * u[0] - u[1]); // dN4/dL1
        dN6[9] = -4. * u[0];                  // dN4/dL2

        dN6[4] = 4. * u[1];                   // dN5/dL1
        dN6[10] = 4. * u[0];                  // dN5/dL2

        dN6[5] = -4. * u[1];                  // dN6/dL1
        dN6[11] = 4. * (1 - u[0] - 2. * u[1]); // dN6/dL2
    }
};

class FemShapeTetra4 : public TemplateShapeFunction<3, 4> 
{
    virtual void computeShapeFunction(const double* x, double* N)
    {
        N[0] = 1. - x[0] - x[1] - x[2];
        N[1] = x[0];
        N[2] = x[1];
        N[3] = x[2];
    }

    virtual void computeGradShapeFunction(const double* /*x*/, double* dNt4)
    {
        //dr
        dNt4[0] = -1.0;
        dNt4[1] = 1.0;
        dNt4[2] = 0.0;
        dNt4[3] = 0.0;

        //ds
        dNt4[4] = -1.0;
        dNt4[5] = 0.0;
        dNt4[6] = 1.0;
        dNt4[7] = 0.0;

        //dt
        dNt4[8] = -1.0;
        dNt4[9] = 0.0;
        dNt4[10] = 0.0;
        dNt4[11] = 1.0;
    }
};

class FemShapeTetra10 : public TemplateShapeFunction<3, 10> 
{
    virtual void computeShapeFunction(const double* x, double* N10)
    {
        N10[0] = 2. * (1 - x[0] - x[1] - x[2]) * (0.5 - x[0] - x[1] - x[2]);
        N10[1] = x[0] * (2. * x[0] - 1);
        N10[2] = x[1] * (2. * x[1] - 1);
        N10[3] = x[2] * (2. * x[2] - 1);
        N10[4] = 4.0 * x[0] * (1.0 - x[0] - x[1] - x[2]);
        N10[5] = 4.0 * x[0] * x[1];
        N10[6] = 4.0 * x[1] * (1.0 - x[0] - x[1] - x[2]);
        N10[7] = 4.0 * x[0] * x[2];
        N10[8] = 4.0 * x[1] * x[2];
        N10[9] = 4.0 * x[2] * (1.0 - x[0] - x[1] - x[2]);
    }

    virtual void computeGradShapeFunction(const double* x, double* dN10)
    {
        dN10[0] = 4.0 * (x[0] + x[1] + x[2]) - 3.0;
        dN10[1] = 4. * x[0] - 1.;
        dN10[2] = 0.0;
        dN10[3] = 0.0;
        dN10[4] = 4.0 * (1.0 - 2.0 * x[0] - x[1] - x[2]);
        dN10[5] = 4.0 * x[1];
        dN10[6] = -4.0 * x[1];
        dN10[7] = 4.0 * x[2];
        dN10[8] = 0.0;
        dN10[9] = -4.0 * x[2];

        dN10[10] =  4. * (x[0] + x[1] + x[2]) - 3.;
        dN10[11] = 0.0;
        dN10[12] = 4. * x[1] - 1.;
        dN10[13] = 0.;
        dN10[14] = -4.0 * x[0];
        dN10[15] = 4.0 * x[0];
        dN10[16] = 4.0 * (1.0 - x[0] - 2.0 * x[1] - x[2]);
        dN10[17] = 0.0;
        dN10[18] = 4.0 * x[2];
        dN10[19] = -4.0 * x[2];

        dN10[20] = 4. * (x[0] + x[1] + x[2]) - 3.;
        dN10[21] = 0.;
        dN10[22] = 0.;
        dN10[23] = 4. * x[2] - 1.;
        dN10[24] = -4.0 * x[0];
        dN10[25] = 0.0;
        dN10[26] = -4.0 * x[1];
        dN10[27] = 4.0 * x[0];
        dN10[28] = 4.0 * x[1];
        dN10[29] = 4.0 * (1.0 - x[0] - x[1] - 2.0 * x[2]);
    }
};

}
