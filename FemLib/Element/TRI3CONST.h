
#pragma once

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "IFemElement.h"

namespace FemLib
{

class TRI3CONST : public IFiniteElement
{
public:
    void assembleNN(MathLib::Matrix<double> &mat) {
        MeshLib::Triangle *ele;
        const size_t n_ele_nodes = 3;
        double a[n_ele_nodes], b[n_ele_nodes], c[n_ele_nodes];
        double nodes_x[n_ele_nodes], nodes_y[n_ele_nodes], nodes_z[n_ele_nodes];
        // xyz
        //for (size_t i=0; i<n_ele_nodes; i++) {
        //    msh->getNodeCoordinates(ele->getNodeID(i), pt);
        //    nodes_x[i] = pt[0];
        //    nodes_y[i] = pt[1];
        //    nodes_z[i] = pt[2];
        //}
        // area
        const double A = 0.5*(nodes_x[0]*(nodes_y[1]-nodes_y[2])+nodes_x[1]*(nodes_y[2]-nodes_y[0])+nodes_x[2]*(nodes_y[0]-nodes_y[1]));
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

        // assemble local EQS
        // Int{w S ph/pt + div(w) K div(p)}dA = Int{w K div(p)}dL
        mat(0,0) = b[0]*b[0] + c[0]*c[0];
        mat(0,1) = b[0]*b[1] + c[0]*c[1];
        mat(0,2) = b[0]*b[2] + c[0]*c[2];
        mat(1,1) = b[1]*b[1] + c[1]*c[1];
        mat(1,2) = b[1]*b[2] + c[1]*c[2];
        mat(2,2) = b[2]*b[2] + c[2]*c[2];
        //local_K *= A;
        // symmetric
        for (size_t i=0; i<n_ele_nodes; i++)
            for (size_t j=0; j<i; j++)
                mat(i,j) = mat(j,i);

    };

    void assembledNdN(MathLib::Matrix<double> &local_K) {
    };
};


}
