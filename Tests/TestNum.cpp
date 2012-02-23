
#include <gtest/gtest.h>

#include <vector>

#include "Base/CodingTools.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"

#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "NumLib/Discrete/DiscreteSystem.h"
#include "NumLib/Discrete/ElementLocalAssembler.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/SparsityBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "NumLib/DomainDecomposition/DomainDecomposition.h"
#include "NumLib/DomainDecomposition/ogs5/par_ddc_group.h"

#include "NumLib/Output/Output.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace NumLib;

TEST(Num, SingleDOF)
{
    DofMapManager dofManagerA;
    size_t dofId1 = dofManagerA.addDoF(10);
    dofManagerA.construct();
    const DofMap *dofMap1 = dofManagerA.getDofMap(dofId1); 

    ASSERT_EQ(dofManagerA.getNumberOfDof(), 1);
    ASSERT_EQ(dofManagerA.getTotalNumberOfDiscretePoints(), 10);
    ASSERT_TRUE(dofMap1!=0);
    ASSERT_EQ(dofMap1->getNumberOfDiscretePoints(), 10);
    ASSERT_EQ(dofMap1->getEqsID(0), 0);
    ASSERT_EQ(dofMap1->getEqsID(9), 9);
};

TEST(Num, NumberingDofByDof)
{
    DofMapManager dofManagerB;
    size_t dofIdB1 = dofManagerB.addDoF(10);
    size_t dofIdB2 = dofManagerB.addDoF(10);
    dofManagerB.construct();
    const DofMap *dofMapB1 = dofManagerB.getDofMap(dofIdB1); 
    const DofMap *dofMapB2 = dofManagerB.getDofMap(dofIdB2); 
    ASSERT_EQ(dofManagerB.getNumberOfDof(), 2);
    ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(dofMapB1->getEqsID(0), 0);
    ASSERT_EQ(dofMapB1->getEqsID(9), 9);
    ASSERT_EQ(dofMapB2->getEqsID(0), 10);
    ASSERT_EQ(dofMapB2->getEqsID(9), 19);
};

TEST(Num, NumberingDofByPoint)
{
    DofMapManager dofManagerB;
    size_t dofIdB1 = dofManagerB.addDoF(10);
    size_t dofIdB2 = dofManagerB.addDoF(10);
    dofManagerB.construct(DofMapManager::BY_POINT);
    const DofMap *dofMapB1 = dofManagerB.getDofMap(dofIdB1); 
    const DofMap *dofMapB2 = dofManagerB.getDofMap(dofIdB2); 
    ASSERT_EQ(dofManagerB.getNumberOfDof(), 2);
    ASSERT_EQ(dofManagerB.getTotalNumberOfDiscretePoints(), 20);
    ASSERT_EQ(dofMapB1->getEqsID(0), 0);
    ASSERT_EQ(dofMapB1->getEqsID(9), 18);
    ASSERT_EQ(dofMapB2->getEqsID(0), 1);
    ASSERT_EQ(dofMapB2->getEqsID(9), 19);
}

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


