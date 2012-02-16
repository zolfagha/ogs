
#include <gtest/gtest.h>

#include <vector>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"


#include "NumLib/Discrete/DiscretizedEQS.h"
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/SparsityBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/DomainDecomposition/DomainDecomposition.h"
#include "NumLib/DomainDecomposition/ogs5/par_ddc_group.h"

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



class TestProblem : public ITransientProblem
{
    //NodeBasedDecomposedDomain* getDomain()
    //{
    //    NodeBasedDecomposedDomain *ddom;

    //    bool isMeshStatic = true;
    //    if (!isMeshStatic) {
    //        MeshLib::IMesh *msh;
    //        DomainDecompostionAlgorithm::decompose(*msh, 1, *ddom);
    //    }
    //    return ddom;
    //};
public:
    int solveTimeStep(const TimeStep &time)
    {
        //NodeBasedDecomposedDomain *ddom = getDomain();
        //ddom->assemble();
        //ddom->solve();
        //ddom->doPost();

        return 0;        
    }
    void initialize()
    {

    }

    double suggestNext(const TimeStep &time_current)
    {
        return .0;
    }
    bool isAwake(const TimeStep &time)
    {
        return true;
    }

};


TEST(Num, test1)
{
    // define problems and solution strategy
    TestProblem problem;

    // pass it to discretization systems
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(problem);

    // start time stepping
    timeStepping.setBeginning(.0);
    timeStepping.solve(100.);
    

    MeshLib::IMixedOrderMesh* msh;
    bool msh_order = false;

    std::vector<ITransientProblem*> problems;
    std::set<std::pair<bool,size_t>> eqs_properties;

    CPARDomain *par;
    CPARDomainGroup doms(*msh, eqs_properties);
    doms.addDomain(par);
    doms.setup();
    doms.solveTimeStep(100.);


}
