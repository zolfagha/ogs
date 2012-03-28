
#pragma once

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/VectorNorms.h"
#include "GeoLib/Shape/Rectangle.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "NumLib/TransientAssembler/ElementLocalAssembler.h"
#include "NumLib/TransientAssembler/TimeEulerElementLocalAssembler.h"
#include "SolutionLib/FemProblem.h"
#include "SolutionLib/SingleStepFEM.h"

#include "Assembler.h"
#include "PorousMedia.h"

using namespace GeoLib;
using namespace MathLib;
using namespace FemLib;
using namespace NumLib;
using namespace MeshLib;
using namespace SolutionLib;
using namespace DiscreteLib;

typedef MathLib::TemplateFunction<GeoLib::Point, double> MyFunction;

class FunctionConcentration : public TemplateTransientMonolithicSystem<2,1>
{
    double _dt;
public:
    FunctionConcentration() {};
    enum Parameters { Velocity=0, Concentration = 1 };
    int solveTimeStep(const TimeStep &ts)
    {
        return 0;
    }
    double suggestNext(const TimeStep &ts) {return ts.getTime()+_dt;};
    bool isAwake(const TimeStep &ts)
    {
        double t = ts.getTime();
        double mod = t - static_cast<int>(t/_dt)*_dt;
        if (mod==.0) return true;
        return false;
    };
    void accept(const TimeStep &time) {};
};
