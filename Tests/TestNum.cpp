
#include <gtest/gtest.h>

#include <vector>

#include "Base/CodingTools.h"
#include "Base/BidirectionalMap.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "DiscreteLib/DiscreteSystem.h"
#include "DiscreteLib/DiscreteLinearEquation.h"
#include "DiscreteLib/DiscreteLinearEquationAssembler.h"
#include "DiscreteLib/ElementLocalAssembler.h"
#include "DiscreteLib/DoF.h"
#include "DiscreteLib/SparsityBuilder.h"
#include "DiscreteLib/ogs5/par_ddc_group.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"

#include "NumLib/Output/Output.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace NumLib;
using namespace DiscreteLib;



TEST(Num, TimeStepping1)
{
    // define problems and solution strategy
    class TestProblem : public ITransientSystem
    {
    private:
        double _t0, _tn, _dt, _v;
    public:
        TestProblem(double t0, double tn, double dt) : _t0(t0), _tn(tn), _dt(dt), _v(.0) {};
        int solveTimeStep(const TimeStep &time) { _v+=1.0; return 0; }
        double suggestNext(const TimeStep &time_current)
        {
            if (time_current.getTime() < _t0) return _t0;
            else return (time_current.getTime()+_dt <= _tn) ? time_current.getTime()+_dt : time_current.getTime();
        }
        bool isAwake(const TimeStep &time)  { return (time.getTime()>=_t0 && time.getTime()<=_tn); }
        void accept(const TimeStep &time) {};
        double getValue() {return _v;};
    };
    TestProblem problem(.0, 100., 10.);

    // pass it to discretization systems
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(problem);
    // start time stepping
    timeStepping.setBeginning(.0);
    size_t n_steps = timeStepping.solve(100.);
    //check
    ASSERT_EQ(10, n_steps);
    ASSERT_EQ(10., problem.getValue());
}

TEST(Num, OGS5DDC)
{
#if 0
    MeshLib::IMixedOrderMesh* msh;
    bool msh_order = false;

    std::vector<ITransientSystem*> problems;
    std::set<std::pair<bool,size_t>> eqs_properties;

    OGS5::CPARDomain *par;
    OGS5::CPARDomainGroup doms(*msh, eqs_properties);
    doms.addDomain(par);
    doms.setup();
    doms.solveTimeStep(100.);
#endif
}

