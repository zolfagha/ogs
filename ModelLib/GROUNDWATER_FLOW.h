
#pragma once

namespace ModeLib
{

typedef int Time;
typedef int ISolution;

class IMathModel
{
public:
    virtual void getFields() = 0;
    virtual void getEquations() = 0; 
};

class Equation
{
    // S ph/pt - div(K grad h) = Q
    // Int{w S ph/pt + grad w dot K grad h}dV = Int{W }
};

class GroundwaterFlowModel : public IMathModel
{
    // unknowns
    //ScalarField head;
    //VectorField velocity;
    // PDE
    Equation pdeVolumeBalance;
    Equation darcyLaw;
    // IC/BC

};


class GROUNDWATER_FLOW
{
public:
    void solve(Time t_n) {

    };

};

}
