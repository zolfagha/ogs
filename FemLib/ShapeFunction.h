
#pragma once

namespace FemLib
{

class IFemShapeFunction
{
public:
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
    virtual size_t getNumerOfNodes() = 0;
    virtual double* getNodeCoordinates(size_t i) = 0;
};

class FemShapeLine2 : public IFemShapeFunction
{
public:
    void computeShapeFunction(const double* pt, double* N)
    {
        N[0] = 1.0 - pt[0];
        N[1] = 1.0 + pt[0];
        for (int i = 0; i < 2; i++)
            N[i] *= 0.5;
    };

    void computeGradShapeFunction(const double* pt, double* dN)
    {
        dN[0] = -0.5;
        dN[1] = 0.5;
    };
};

class FemShapeLine3 : public IFemShapeFunction
{
public:
    void computeShapeFunction(const double* pt, double* N)
    {
        N[0] = 0.5 * pt[0] * (pt[0] - 1.0);
        N[1] = 0.5 * pt[0] * (pt[0] + 1.0);
        N[2] = 1.0 - pt[0] * pt[0];
    };

    void computeGradShapeFunction(const double* pt, double* dN)
    {
        dN[0] = pt[0] - 0.5;
        dN[1] = pt[0] + 0.5;
        dN[2] = -2.0 * pt[0];
    };
};

class FemShapeQuad4 : public IFemShapeFunction
{
public:
    void computeShapeFunction(const double* pt, double* N)
    {
        N[0] = (1.0 + pt[0]) * (1.0 + pt[1]);
        N[1] = (1.0 - pt[0]) * (1.0 + pt[1]);
        N[2] = (1.0 - pt[0]) * (1.0 - pt[1]);
        N[3] = (1.0 + pt[0]) * (1.0 - pt[1]);
        for (int i = 0; i < 4; i++)
            N[i] *= 0.25;
    };

    void computeGradShapeFunction(const double* u, double* dN4)
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

class FemShapeQuad8 : public IFemShapeFunction
{
public:
    void computeShapeFunction(const double* pt, double* N)
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

    void computeGradShapeFunction(const double* pt, double* dN8)
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

class FemShapeQuad9 : public IFemShapeFunction
{
public:
    void computeShapeFunction(const double* pt, double* N)
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

    void computeGradShapeFunction(const double* u, double* dN9)
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

class FemShapeHex8 : public IFemShapeFunction
{
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
};

class FemShapeHex20 : public IFemShapeFunction
{
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
};


class FemShapeTriangle3 : public IFemShapeFunction
{
    virtual void computeShapeFunction(double* u, double* N)
    {
        N[0] = 1. - u[0] - u[1];
        N[1] = u[0];
        N[2] = u[1];
    };
    virtual void computeGradShapeFunction(double* pt, double* dN3)
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

class FemShapeTriangle6 : public IFemShapeFunction
{
    virtual void computeShapeFunction(double* u, double* N6)
    {
        N6[0] = 2. * (1. - u[0] - u[1]) * (0.5 - u[0] - u[1]);
        N6[1] = u[0] * (2. * u[0] - 1.);
        N6[2] = u[1] * (2. * u[1] - 1.);
        N6[3] = 4. * u[0] * (1. - u[0] - u[1]);
        N6[4] = 4. * u[0] * u[1];
        N6[5] = 4. * u[1] * (1. - u[0] - u[1]);
    };
    virtual void computeGradShapeFunction(double* u, double* dN6)
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

class FemShapeTetra4 : public IFemShapeFunction
{
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
};

class FemShapeTetra10 : public IFemShapeFunction
{
    virtual void computeShapeFunction(double* pt, double* N) = 0;
    virtual void computeGradShapeFunction(double* pt, double* dN) = 0;
};



};
