
#include "TRI3CONST.h"


#include "MathLib/LinAlg/Dense/Matrix.h"
#include "GeoLib/Core/Point.h"
#include "IFemElement.h"

namespace FemLib
{

void TRI3CONST::configure( MeshLib::IMesh * msh, MeshLib::Triangle * e )
{
    _msh = msh;
    _ele = e;
    double nodes_x[3], nodes_y[3];
    // xyz
    for (size_t i=0; i<3; i++) {
        const GeoLib::Point *pt = _msh->getNode(_ele->getNodeID(i))->getData();
        nodes_x[i] = pt->getData()[0];
        nodes_y[i] = pt->getData()[1];
    }
    // area
    A = _ele->getArea();
    // set a,b,c
    a[0] = 0.5/A*(nodes_x[1]*nodes_y[2]-nodes_x[2]*nodes_y[1]);
    b[0] = 0.5/A*(nodes_y[1]-nodes_y[2]);
    c[0] = 0.5/A*(nodes_x[2]-nodes_x[1]);
    a[1] = 0.5/A*(nodes_x[2]*nodes_y[0]-nodes_x[0]*nodes_y[2]);
    b[1] = 0.5/A*(nodes_y[2]-nodes_y[0]);
    c[1] = 0.5/A*(nodes_x[0]-nodes_x[2]);
    a[2] = 0.5/A*(nodes_x[0]*nodes_y[1]-nodes_x[1]*nodes_y[0]);
    b[2] = 0.5/A*(nodes_y[0]-nodes_y[1]);
    c[2] = 0.5/A*(nodes_x[1]-nodes_x[0]);
}

double TRI3CONST::interpolate(double *x, double *nodal_values)
{
    double v = 0;
    for (size_t i=0; i<3; i++)
        v += (a[i]+b[i]*x[0]+c[i]*x[1])*nodal_values[i];
    return v;
}

/// compute an matrix M = Int{W^T F N} dV
void TRI3CONST::computeIntTestShape( Fscalar f, MathLib::Matrix<double> &mat)
{
    const double v = f(0);
    mat(0,0) = 1.0;
    mat(0,1) = 0.5;
    mat(0,2) = 0.5;
    mat(1,1) = 1.0;
    mat(1,2) = 0.5;
    mat(2,2) = 1.0;
    mat *= v*A/6.0;
    // make symmetric
    for (size_t i=0; i<3; i++)
        for (size_t j=0; j<i; j++)
            mat(i,j) = mat(j,i);
}

/// compute an matrix M = Int{W^T F dN} dV
void TRI3CONST::computeIntTestDShape( Fvector f, MathLib::Matrix<double> &mat)
{
    double *v = f(0);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            mat(i,j) = v[0]*b[j] + v[1]*c[j];
    mat *= A/3.0;
}

/// compute an matrix M = Int{dW^T F dN} dV
void TRI3CONST::computeIntDTestDShape( Fscalar f, MathLib::Matrix<double> &mat)
{
    const double v = f(0);
    mat(0,0) = b[0]*b[0] + c[0]*c[0];
    mat(0,1) = b[0]*b[1] + c[0]*c[1];
    mat(0,2) = b[0]*b[2] + c[0]*c[2];
    mat(1,1) = b[1]*b[1] + c[1]*c[1];
    mat(1,2) = b[1]*b[2] + c[1]*c[2];
    mat(2,2) = b[2]*b[2] + c[2]*c[2];
    mat *= v*A;
    // make symmetric
    for (size_t i=0; i<3; i++)
        for (size_t j=0; j<i; j++)
            mat(i,j) = mat(j,i);
}

}
