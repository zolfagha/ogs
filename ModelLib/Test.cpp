
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "GROUNDWATER_FLOW.h"

using namespace ModeLib;

class TimeStep {
public:
    bool next();
    Time getTime();
};

template <typename T>
class Field
{
public:
    virtual T& get(double*pt) = 0;
};

typedef Field<double> ScalarField;
typedef Field<double*> VectorField;
typedef Field<MathLib::Matrix<double> > TensorField;


struct PorousMediumProperty
{
    MathLib::Matrix<double> K;
    double n;
};

struct FluidProperty
{
    double rho;
    double cp;
    MathLib::Matrix<double> lambda;
};


int main(int argc, char* argv[]) 
{

    ScalarField *head; //fem, nodal
    VectorField *velocity; //fem, integration pt
    TensorField *K; // homogeneous, element hetero, function(x), function(h, T)
    Field<PorousMediumProperty> *poro; // function()
    Field<FluidProperty> *fluid; // function()

    GROUNDWATER_FLOW gw; //gw(poro)

    TimeStep time;

    while (time.next()) {
        Time t = time.getTime();
        gw.solve(t);

        //while (!converge(head, T)) {
        // materials.update(head, velocity, T);
        // head = gw.solveHead(t, materials, head_n);
        // velocity = FluidMomentum.solveVelocity(head);
        // T = Heat.solveT(t, T_n, velocity);
        //}

        //while (!converge) {
        // {head, T} = GwT.solve(t, materials, velocity);
        // vel = FluidMomentum.solve(head);
        //}

        // disp = df.solve();

        //output
    }

    return 0;
};

