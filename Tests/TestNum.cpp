
#include <gtest/gtest.h>

#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "NumLib/Discrete/DiscretizedEQS.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/SparsityBuilder.h"

#include "NumLib/Coupling/Clock.h"
#include "NumLib/Coupling/TransientSystems.h"
#include "NumLib/Coupling/AsyncPartSolution.h"
#include "NumLib/Coupling/PartitionedAlgorithm.h"
#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"

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

TEST(Num, testDis1)
{    
    //mesh
    IMesh* msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    MeshLib::TopologyNode2Nodes topo_node2nodes(msh);
    //define dof
    DofMapManager dofManager;
    size_t dofId = dofManager.addDoF(msh->getNumberOfNodes());
    dofManager.construct();
    //sparse table
    RowMajorSparsity sparse;
    createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparse);

    // construct discrete eqs
    SparseLinearEquations eqs;
    eqs.create(sparse.size(), &sparse);
    eqs.getOption().solver_type = SparseLinearEquations::SolverCG;
    eqs.getOption().precon_type = SparseLinearEquations::NONE;

    //
    const DofMap *dofMap = dofManager.getDofMap(dofId);
    for (size_t i=0; i<msh->getNumberOfElements(); i++) {
        const IElement* e = msh->getElemenet(i);
        std::vector<size_t> ele_node_ids, local_dofmap;
        e->getNodeIDList(ele_node_ids);
        dofMap->getListOfEqsID(ele_node_ids, local_dofmap);
        MathLib::Matrix<double> localK(local_dofmap.size(),local_dofmap.size());
        eqs.addA(local_dofmap, localK);
    }

    eqs.solve();
}


// 2a + 2b + .3c = 6.9
// 3a + 5b + .2c = 13.6
// .5a + .3b + 3c = 10.1
// A. a=1, b=2, c=3
class WeakCouplingEQS1 : public TemplateMonolithicProblem<2,1>
{
public:
    enum Parameters { a=2, b = 0, c = 1 };
    int solve()
    {
        double vb = _vec_parameters[WeakCouplingEQS1::b]->eval(.0);
        double vc = _vec_parameters[WeakCouplingEQS1::c]->eval(.0);
        double va = 1./2.*(6.9 - 2.*vb - 0.3*vc);
        if (_vec_parameters[WeakCouplingEQS1::a]!=0)
            delete _vec_parameters[WeakCouplingEQS1::a];
        _vec_parameters[WeakCouplingEQS1::a] = new MathLib::FunctionConstant<double,double>(va);
        return 0;
    }
};

class WeakCouplingEQS2 :  public TemplateMonolithicProblem<2,1>
{
public:
    enum Parameters { a = 0, b = 2, c = 1 };
    int solve()
    {
        double va = _vec_parameters[WeakCouplingEQS2::a]->eval(.0);
        double vc = _vec_parameters[WeakCouplingEQS2::c]->eval(.0);
        double vb = 1./5.*(13.6-3*va-0.2*vc);
        if (_vec_parameters[WeakCouplingEQS2::b]!=0)
            delete _vec_parameters[WeakCouplingEQS2::b];
        _vec_parameters[WeakCouplingEQS2::b] = new MathLib::FunctionConstant<double,double>(vb);
        return 0;
    }
};

class WeakCouplingEQS3 : public TemplateMonolithicProblem<2,1>
{
public:
    enum Parameters { a = 0, b = 1, c = 2 };
    int solve()
    {
        double va = _vec_parameters[a]->eval(.0);
        double vb = _vec_parameters[b]->eval(.0);
        double vc = 1./3.*(10.1-0.5*va-0.3*vb);
        if (_vec_parameters[c]!=0)
            delete _vec_parameters[c];
        _vec_parameters[c] = new MathLib::FunctionConstant<double,double>(vc);
        return 0;
    }
private:
    double va, vb, vc;
};

#if 0
// 2a + 2b + 3c = 15
// 3a + 5b + 2c = 19
// 5a + 3b + 3c = 20
// A. a=1, b=2, c=3
class StrongCouplingEQS1 : public IMonolithicProblem
{
public:
    enum InputParameters { b = 0, c = 1 };
    enum OutputParameters { a = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS1() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double vb = _vec_in_var[EQS1::b]->eval(.0);
        double vc = _vec_in_var[EQS1::c]->eval(.0);
        double va = 1./2.*(15. - 2.*vb - 3.*vc);
        if (_vec_out_var[EQS1::a]!=0)
            delete _vec_out_var[EQS1::a];
        _vec_out_var[EQS1::a] = new MathLib::FunctionConstant<double,double>(va);
        return 0;
    }
};

class StrongCouplingEQS2 : public IMonolithicProblem
{
public:
    enum InputParameters { a = 0, c = 1 };
    enum OutputParameters { b = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS2() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double va = _vec_in_var[EQS2::a]->eval(.0);
        double vc = _vec_in_var[EQS2::c]->eval(.0);
        double vb = 1./5.*(19.-3*va-2.*vc);
        if (_vec_out_var[EQS2::b]!=0)
            delete _vec_out_var[EQS2::b];
        _vec_out_var[EQS2::b] = new MathLib::FunctionConstant<double,double>(vb);
        return 0;
    }
};

class StrongCouplingEQS3 : public IMonolithicProblem
{
public:
    enum InputParameters { a = 0, b = 1 };
    enum OutputParameters { c = 0 };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    StrongCouplingEQS3() {
        _vec_in_var.resize(getNumberOfInputVarameters(), 0);
        _vec_out_var.resize(getNumberOfOutputParameters(), 0);
    }
    int solve()
    {
        double va = _vec_in_var[a]->eval(.0);
        double vb = _vec_in_var[b]->eval(.0);
        double vc = 1./3.*(20.-5*va-3.*vb);
        if (_vec_out_var[c]!=0)
            delete _vec_out_var[c];
        _vec_out_var[c] = new MathLib::FunctionConstant<double,double>(vc);
        return 0;
    }
private:
    double va, vb, vc;
};
#endif


TEST(Coupling, SteadyCouplingCheck1)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    eqs1.setParameter(WeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(WeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(WeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    {
        // correct
        PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
        part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
        part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
        part1.addParameter("c");
        part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
        part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
        part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
        part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

        PartitionedProblem part2(BlockJacobiMethod(1.e-4, 100));
        part2.addParameter("a", part1, part1.getParameterID("a"));
        part2.addParameter("b", part1, part1.getParameterID("b"));
        part2.addParameter("c", eqs3, WeakCouplingEQS3::c);
        part2.connectInput("a", eqs3, WeakCouplingEQS3::a);
        part2.connectInput("b", eqs3, WeakCouplingEQS3::b);
        part2.connectInput("c", part1, part1.getParameterID("c"));

        ASSERT_TRUE(part2.check());
    }

    {
        // no source defined for c
        PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
        part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
        part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
        part1.addParameter("c");
        part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
        part1.connectInput("a", eqs3, WeakCouplingEQS3::a);
        part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
        part1.connectInput("b", eqs3, WeakCouplingEQS3::b);
        part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
        part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

        ASSERT_TRUE(part1.check());
    }

    {
        // missing connections eqs2:c
        PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
        part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
        part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
        part1.addParameter("c", eqs3, WeakCouplingEQS3::c);
        part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
        part1.connectInput("a", eqs3, WeakCouplingEQS3::a);
        part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
        part1.connectInput("b", eqs3, WeakCouplingEQS3::b);
        part1.connectInput("c", eqs1, WeakCouplingEQS1::c);

        ASSERT_FALSE(part1.check());
    }

}


TEST(Coupling, SteadyCouplingJacobi)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    eqs1.setParameter(WeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(WeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(WeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    PartitionedProblem part1(BlockJacobiMethod(1.e-4, 100));
    part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
    part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
    part1.addParameter("c");
    part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
    part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
    part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
    part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

    PartitionedProblem part2(BlockJacobiMethod(1.e-4, 100));
    part2.addParameter("a", part1, part1.getParameterID("a"));
    part2.addParameter("b", part1, part1.getParameterID("b"));
    part2.addParameter("c", eqs3, WeakCouplingEQS3::c);
    part2.connectInput("a", eqs3, WeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, WeakCouplingEQS3::b);
    part2.connectInput("c", part1, part1.getParameterID("c"));

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    ASSERT_NEAR(1., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(2., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(3., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

TEST(Coupling, SteadyCouplingSeidel)
{
    WeakCouplingEQS1 eqs1;
    WeakCouplingEQS2 eqs2;
    WeakCouplingEQS3 eqs3;

    eqs1.setParameter(WeakCouplingEQS1::a, new MathLib::FunctionConstant<double,double>(.0));
    eqs2.setParameter(WeakCouplingEQS2::b, new MathLib::FunctionConstant<double,double>(.0));
    eqs3.setParameter(WeakCouplingEQS3::c, new MathLib::FunctionConstant<double,double>(.0));

    PartitionedProblem part1(BlockGaussSeidelMethod(1.e-5, 100));
    part1.addParameter("a", eqs1, WeakCouplingEQS1::a);
    part1.addParameter("b", eqs2, WeakCouplingEQS2::b);
    part1.addParameter("c");
    part1.connectInput("b", eqs1, WeakCouplingEQS1::b);
    part1.connectInput("c", eqs1, WeakCouplingEQS1::c);
    part1.connectInput("a", eqs2, WeakCouplingEQS2::a);
    part1.connectInput("c", eqs2, WeakCouplingEQS2::c);

    PartitionedProblem part2(BlockGaussSeidelMethod(1.e-5, 100));
    part2.addParameter("a", part1, part1.getParameterID("a"));
    part2.addParameter("b", part1, part1.getParameterID("b"));
    part2.addParameter("c", eqs3, WeakCouplingEQS3::c);
    part2.connectInput("a", eqs3, WeakCouplingEQS3::a);
    part2.connectInput("b", eqs3, WeakCouplingEQS3::b);
    part2.connectInput("c", part1, part1.getParameterID("c"));

    ASSERT_TRUE(part2.check());
    part2.solve();

    const double epsilon = 1.e-3;
    ASSERT_NEAR(1., part2.getParameter(part2.getParameterID("a"))->eval(0), epsilon);
    ASSERT_NEAR(2., part2.getParameter(part2.getParameterID("b"))->eval(0), epsilon);
    ASSERT_NEAR(3., part2.getParameter(part2.getParameterID("c"))->eval(0), epsilon);
}

//TEST(Coupling, TransientCoupling)
//{
//    std::vector<ITransientSystem*> vec_systems;
//
//    Clock clock;
//    for (size_t i=0; i<vec_systems.size(); i++)
//        clock.addTransientSystem(vec_systems[i]);
//
//    TimeStep t0=.0, t_end=1.0;
//    clock.setBeginning(t0);
//
//    // start clock
//    clock.moveForwardUntill(t_end); 
//}

