
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
#include "NumLib/Discrete/DoF.h"
#include "NumLib/Discrete/SparsityBuilder.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"
#include "NumLib/DomainDecomposition/DomainDecomposition.h"
#include "NumLib/DomainDecomposition/ogs5/par_ddc_group.h"

#include "NumLib/Solution/ISolution.h"
#include "NumLib/Solution/IProblem.h"

#include "NumLib/Output/Output.h"

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
    MeshLib::IMixedOrderMesh* msh;
    bool msh_order = false;

    std::vector<ITransientSystem*> problems;
    std::set<std::pair<bool,size_t>> eqs_properties;

    OGS5::CPARDomain *par;
    OGS5::CPARDomainGroup doms(*msh, eqs_properties);
    doms.addDomain(par);
    doms.setup();
    doms.solveTimeStep(100.);
}


class GWAssembly: public AbstractTimeODEFemElementAssembler
{
private:
    MathLib::IFunction<double, double*>* _matK;
public:
    GWAssembly(MathLib::IFunction<double, double*> &mat) : _matK(&mat) {};
    //explicit GWAssembly(FemLib::FemNodalFunctionScalar &fem, MathLib::IFunction<double, double*> &mat) : AbstractTimeODEFemElementAssembler(fem), _matK(&mat) {};

protected:
    void assemblyFE(const TimeStep &time, FemLib::IFiniteElement &fe, MathLib::DenseLinearEquations::MatrixType &localM, MathLib::DenseLinearEquations::MatrixType &localK,  MathLib::DenseLinearEquations::VectorType &localF)
    {
        localM = .0;
        localF.resize(localF.size(), .0);
        fe.integrateDWxDN(_matK, localK);
        localK *= -1;
    }
};

class GWFemTestSystem : public ITransientSystem
{
public:

    GWFemTestSystem()
    {
        Base::zeroObject(rec, K, head, _discrete);
    }
    virtual ~GWFemTestSystem()
    {
        Base::releaseObject(rec, K, head);
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
        head = _solHead->getCurrentSolution(0);
        return 0;
    }


    //#Define a problem
    void define()
    {
        //geometry
        rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = rec->getLeft();
        Polyline* poly_right = rec->getRight();
        //mesh
        MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
        //mat
        K = new MathLib::FunctionConstant<double, double*>(1.e-11);
        //equations
        GWAssembly ele_eqs(*K) ;
        //IVBV problem
        _problem = new FemIVBVProblem<GWAssembly>(*msh, ele_eqs);
        size_t headId = _problem->createField(PolynomialOrder::Linear);
        _problem->addDirichletBC(headId, *poly_right, MathLib::FunctionConstant<double, GeoLib::Point>(.0));
        _problem->addNeumannBC(headId, *poly_left, MathLib::FunctionConstant<double, GeoLib::Point>(-1e-5));
        _problem->setTimeSteppingFunction(TimeStepFunctionConstant(.0, 100.0, .0));
        //solution algorithm
        _solHead = new TimeEulerSpaceFemLinearAlgorithm<FemIVBVProblem<GWAssembly>, GWAssembly>(*_problem);


        head = _problem->getField(headId);
        //vel = new FEMIntegrationPointFunctionVector2d(msh);
    }


    void setDiscreteSystem(DiscreteSystem& dis)
    {
        _discrete = &dis;
    }


private:
    FemIVBVProblem<GWAssembly>* _problem;
    TimeEulerSpaceFemLinearAlgorithm<FemIVBVProblem<GWAssembly>, GWAssembly>* _solHead; 
    Rectangle *rec;

    MathLib::IFunction<double, double*> *K;
    FemNodalFunctionScalar *head;
    DiscreteSystem* _discrete;
};





TEST(Num, Discrete1)
{
    {
    // define problems
    GWFemTestSystem gwProblem;
    gwProblem.define();
    //// create discrete systems according to mesh
    //MeshLib::IMesh *msh = gwProblem.getMesh();
    //DiscreteSystem dis(*msh, *new MathLib::SparseLinearEquations());
    //// 
    //gwProblem.setDiscreteSystem(dis);

    // start time stepping
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(gwProblem);

    timeStepping.setBeginning(.0);
    timeStepping.solve(100.);
    }


}

