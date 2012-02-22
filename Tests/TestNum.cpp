
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

#include "NumLib/Solution/ISolution.h"
#include "NumLib/Solution/IProblem.h"

#include "NumLib/Output/Output.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace FemLib;
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
    MeshLib::TopologyNode2NodesConnectedByEdges topo_node2nodes(msh);
    //define dof
    DofMapManager dofManager;
    size_t dofId = dofManager.addDoF(msh->getNumberOfNodes());
    dofManager.construct();
    //sparse table
    RowMajorSparsity sparse;
    SparsityBuilder::createRowMajorSparsityFromNodeConnectivity(topo_node2nodes, sparse);

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




class TestProblem : public ITransientSystem
{
public:
    int solveTimeStep(const TimeStep &time)
    {

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

    void accept(const TimeStep &time) {};
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


class GWAssembler: public AbstractTimeODEFemElementAssembler
{
private:
    MathLib::IFunction<double, double*>* _matK;
    LagrangianFeObjectContainer* _feObjects;
public:
    GWAssembler(LagrangianFeObjectContainer &feObjects, MathLib::IFunction<double, double*> &mat) : _feObjects(&feObjects), _matK(&mat) 
    {
    };

//protected:
    void assembly(const TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::MatrixType &localM, MathLib::DenseLinearEquations::MatrixType &localK,  MathLib::DenseLinearEquations::VectorType &localF)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);

        //localM = .0;
        //localF.resize(localF.size(), .0);
        fe->integrateDWxDN(_matK, localK);
    }
};

class GWFemTestSystem : public ITransientSystem
{
    typedef FemIVBVProblem<GWAssembler> GWFemProblem;
    typedef SingleStepLinearFEM<TimeEulerElementAssembler, MathLib::SparseLinearEquations, FemIVBVProblem<GWAssembler>> SolutionForHead;
public:
    GWFemTestSystem()
    {
        Base::zeroObject(_rec, _head, _problem);
    }

    virtual ~GWFemTestSystem()
    {
        Base::releaseObject(_rec, _head, _problem);
    }

    double suggestNext(const TimeStep &time_current)
    {
        return _solHead->suggestNext(time_current);
    }
    bool isAwake(const TimeStep &time)
    {
        return _solHead->isAwake(time);
    }

    int solveTimeStep(const TimeStep &time)
    {
        _solHead->solveTimeStep(time);
        _head = _solHead->getCurrentSolution(0);
        return 0;
    }

    void accept(const TimeStep &time) 
    {
        _solHead->accept(time);
    };

    //#Define a problem
    void define(DiscreteSystem &dis, MathLib::IFunction<double, double*> &K)
    {
        MeshLib::IMesh *msh = dis.getMesh();
        size_t nnodes = msh->getNumberOfNodes();
        _feObjects = new LagrangianFeObjectContainer(*msh);
        //equations
        GWAssembler ele_eqs(*_feObjects, K) ;
        //IVBV problem
        _problem = new GWFemProblem(*dis.getMesh(), GWAssembler(ele_eqs));
        //BC
        size_t headId = _problem->createField(PolynomialOrder::Linear);
        _head = _problem->getField(headId);
        _rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = _rec->getLeft();
        Polyline* poly_right = _rec->getRight();
        _problem->setIC(headId, *_head);
        _problem->addDirichletBC(headId, *poly_right, false, MathLib::FunctionConstant<double, GeoLib::Point>(.0));
        _problem->addNeumannBC(headId, *poly_left, false, MathLib::FunctionConstant<double, GeoLib::Point>(-1e-5));
        //transient
        _problem->setTimeSteppingFunction(TimeStepFunctionConstant(.0, 100.0, 10.0));
        //solution algorithm
        _solHead = new SolutionForHead(dis, *_problem);
        _solHead->getTimeODEAssembler()->setTheta(1.0);
        MathLib::SparseLinearEquations* linear_solver = _solHead->getLinearEquationSolver();
        linear_solver->getOption().solver_type = SparseLinearEquations::SolverCG;
        linear_solver->getOption().precon_type = SparseLinearEquations::NONE;


        //vel = new FEMIntegrationPointFunctionVector2d(msh);
    }

    FemNodalFunctionScalar* getCurrentHead()
    {
        return _head;
    }

private:
    GWFemProblem* _problem;
    SolutionForHead* _solHead; 
    Rectangle *_rec;
    FemNodalFunctionScalar *_head;
    LagrangianFeObjectContainer* _feObjects;

    DISALLOW_COPY_AND_ASSIGN(GWFemTestSystem);
};





TEST(Num, DiscreteSingleStep1)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::FunctionConstant<double, double*> K(1.e-11);
    // BC
    // define problems
    GWFemTestSystem gwProblem;
    gwProblem.define(dis, K);

    gwProblem.solveTimeStep(TimeStep(1.0, 1.0));

    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h->getNodalValues(), h->getNumberOfNodes());

    Base::releaseObject(msh);
}

TEST(Num, DiscreteTimeStep1)
{
    // create a discrete system
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
    DiscreteSystem dis(*msh);
    // mat
    MathLib::FunctionConstant<double, double*> K(1.e-11);
    // BC
    // define problems
    GWFemTestSystem gwProblem;
    gwProblem.define(dis, K);

    // start time stepping
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(gwProblem);
    timeStepping.setBeginning(.0);
    timeStepping.solve(100.);


    std::vector<double> expected(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) expected[i] = 2.e+6;
        if (i%3==1) expected[i] = 1.e+6;
        if (i%3==2) expected[i] = 0.e+6;
    }

    FemNodalFunctionScalar* h = gwProblem.getCurrentHead();

    ASSERT_DOUBLE_ARRAY_EQ(&expected[0], h->getNodalValues(), h->getNumberOfNodes());


    Base::releaseObject(msh);
}
