
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
#include "NumLib/Coupling/CouplingSolution.h"
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


class BiotEQS : public IMonolithicProblem
{
public:
    enum InputParameter
    {
        p = 0,
        T = 1
    };

    enum OutputParameter
    {
        u = 0
    };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};

    BiotEQS() {
        _vec_in_var.resize(getNumberOfInputVarameters());
        _vec_out_var.resize(getNumberOfOutputParameters());
    }

    int solve()
    {
        // update u
        return 0;
    }
};

class LiquidEQS : public IMonolithicProblem
{
public:
    enum InputParameters
    {
        u = 0,
        T = 1
    };
    enum OutputParameters
    {
        p = 0,
        v = 1
    };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 2;};
    LiquidEQS() {
        _vec_in_var.resize(getNumberOfInputVarameters());
        _vec_out_var.resize(getNumberOfOutputParameters());
    }
    int solve()
    {
        // update u
        return 0;
    }

};

class HeatEQS : public IMonolithicProblem
{
public:
    enum InputParameters
    {
        v = 0,
        p = 1
    };
    enum OutputParameters
    {
        T = 0
    };
    size_t getNumberOfInputVarameters() const {return 2;};
    size_t getNumberOfOutputParameters() const {return 1;};
    HeatEQS() {
        _vec_in_var.resize(getNumberOfInputVarameters());
        _vec_out_var.resize(getNumberOfOutputParameters());
    }
    int solve()
    {
        // update u
        return 0;
    }
};


TEST(Coupling, SteadyCoupling1)
{
    BiotEQS biotEQS;
    LiquidEQS liquid;
    HeatEQS heat;

    biotEQS.setInitial(BiotEQS::u, new MathLib::FunctionConstant<double,double>(.0));
    liquid.setInitial(LiquidEQS::p, new MathLib::FunctionConstant<double,double>(1.0e+6));
    liquid.setInitial(LiquidEQS::v, new MathLib::FunctionConstant<double,double>(.0));
    heat.setInitial(HeatEQS::T, new MathLib::FunctionConstant<double,double>(25.0));

    PartitionedProblem partPU(&BlockJacobiMethod());
    partPU.add("u", &biotEQS, BiotEQS::u);
    partPU.add("p", &liquid, LiquidEQS::p);
    partPU.add("T");
    partPU.connectInput("u", &liquid, LiquidEQS::u);
    partPU.connectInput("p", &biotEQS, BiotEQS::p);
    partPU.connectInput("T", &liquid, LiquidEQS::T);
    partPU.connectInput("T", &biotEQS, BiotEQS::T);

    PartitionedProblem partTV(&BlockJacobiMethod());
    partTV.add("v", &partPU, partPU.getVariableID("v"));
    partTV.add("p", &partPU, partPU.getVariableID("p"));
    partTV.add("T", &heat, HeatEQS::T);
    partTV.connectInput("v", &heat, HeatEQS::v);
    partTV.connectInput("p", &heat, HeatEQS::p);
    partTV.connectInput("T", &partPU, partPU.getVariableID("T"));



    // v = P()
    // T = Heat(v)

    ASSERT_FALSE(partTV.check());
}


#if 0
TEST(Coupling, SteadyCoupling2)
{
    SharedVariables vars;
    PartitionedProblem partSol1(&BlockJacobiMethod());
    PartitionedProblem partSol2(&BlockJacobiMethod());
    MonolithicProblem mono1;
    MonolithicProblem mono2;
    MonolithicProblem mono3;
    partSol1.addChildren(&mono1);
    partSol1.addChildren(&partSol2);
    partSol2.addChildren(&mono2);
    partSol2.addChildren(&mono3);
    partSol1.solve(&vars);
}
#endif

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

