
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
#include "NumLib/DomainDecomposition/DomainDecomposition.h"
#include "NumLib/DomainDecomposition/ogs5/par_ddc_group.h"

#include "NumLib/Solution/ISolution.h"
#include "NumLib/Solution/TransientFEModel.h"

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




class TestProblem : public ITransientProblem
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

    std::vector<ITransientProblem*> problems;
    std::set<std::pair<bool,size_t>> eqs_properties;

    OGS5::CPARDomain *par;
    OGS5::CPARDomainGroup doms(*msh, eqs_properties);
    doms.addDomain(par);
    doms.setup();
    doms.solveTimeStep(100.);
}


class GWFemTestProblem : public ITransientProblem
{
public:
    GWFemTestProblem()
    {
        Base::zeroObject(rec, msh, K, head, vel, _discrete);
    }
    virtual ~GWFemTestProblem()
    {
        Base::releaseObject(rec, msh, K, head, vel);
        Base::releaseObjectsInStdVector(vec_bc1);
        Base::releaseObjectsInStdVector(vec_bc2);
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

    //#Define a problem
    void define()
    {
        //geometry
        rec = new Rectangle(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
        Polyline* poly_left = rec->getLeft();
        Polyline* poly_right = rec->getRight();
        //mesh
        this->msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);
        //discretization
        head = new FemNodalFunctionScalar(msh, PolynomialOrder::Linear);
        head_n = new FemNodalFunctionScalar(msh, PolynomialOrder::Linear);
        vel = new FEMIntegrationPointFunctionVector2d(msh);
        //bc
        vec_bc1.push_back(new FemDirichletBC<double>(head, poly_right, new MathLib::FunctionConstant<double, GeoLib::Point>(.0), new DiagonalizeMethod())); //TODO should BC objects be created by fe functions?
        vec_bc2.push_back(new FemNeumannBC<double, double>(head, poly_left, new MathLib::FunctionConstant<double, GeoLib::Point>(-1e-5)));

        K = new MathLib::FunctionConstant<double, double*>(1.e-11);
    }

    int solveTimeStep(const TimeStep &time)
    {
        class GWAssembly: public ITimeODEFemElementAssembler
        {
        public:
            void assembly(TimeStep &time, FemLib::IFiniteElement &fe, MathLib::DenseLinearEquations::MatrixType &localM, MathLib::DenseLinearEquations::MatrixType &localK,  MathLib::DenseLinearEquations::VectorType &localF)
            {
                MathLib::FunctionConstant<double, double*> K(1.e-11); //where do you get this material properties?

                localM = .0;
                localF.resize(localF.size(), .0);
                fe.integrateDWxDN(&K, localK);
                localK *= -1;
            }
        };


        TimeEulerSpFEMLinearSolution<GWAssembly> sol; 
        head = sol.solve(time, *head_n);

        return 0;
    }

    void setDiscreteSystem(DiscreteSystem& dis)
    {
        _discrete = &dis;
    }

    IMixedOrderMesh* getMesh() const {return msh;};

private:
    Rectangle *rec;
    IMixedOrderMesh *msh;
    MathLib::IFunction<double, double*> *K;
    FemNodalFunctionScalar *head;
    FemNodalFunctionScalar *head_n;
    FEMIntegrationPointFunctionVector2d *vel;
    std::vector<FemDirichletBC<double>*> vec_bc1;
    std::vector<FemNeumannBC<double, double>*> vec_bc2;
    DiscreteSystem* _discrete;
};

enum TimeUnit
{
    Second,
    Minute,
    Day,
    Week,
    Month,
    Year
};

class ITimeStepFunction
{

};


class TimeStepFunctionConstant : public ITimeStepFunction
{
public:
    TimeStepFunctionConstant(double dt) {};
};

enum TimeType
{
    ALL_STEPS
};

class IOutput
{
public:
    IOutput(TimeType time_tpye) {};

    void addNodal(int);
    void addElemental(int);
};

class PVD : public IOutput
{
public:
    PVD(TimeType time_tpye) : IOutput(time_tpye) {};
};


class TECPLOT : public IOutput
{
public:
    TECPLOT(TimeType time_tpye) : IOutput(time_tpye) {};
    void setGeometry(GeoLib::GeoObject& geo);
};

class TransientFemModel
{
public:
    void setMesh(MeshLib::IMesh &msh);
    void addBC(int var_type, GeoLib::GeoObject& geo, MathLib::IFunction<double, double>& f);
    void addST(int var_type, GeoLib::GeoObject& geo, MathLib::IFunction<double, double>& f);
    void setTransient(TimeUnit unit, double t0, double tn, ITimeStepFunction &f);
    void setLinearSolver(std::string& solver, double epsilon, size_t max_itr);
    void addOutput(IOutput& out);
    void construct();
    void doAnalysis();
    MathLib::IFunction<double, double>* getResult(int);
};

class GROUNDWATER_FLOW : public TransientFemModel
{
private:
    ISolutionAlgorithm* _solution;

    class GWAssembly: public ITimeODEFemElementAssembler
    {
    public:
        void assembly(TimeStep &time, FemLib::IFiniteElement &fe, MathLib::DenseLinearEquations::MatrixType &localM, MathLib::DenseLinearEquations::MatrixType &localK,  MathLib::DenseLinearEquations::VectorType &localF)
        {
            MathLib::FunctionConstant<double, double*> K(1.e-11); //where do you get this material properties?

            localM = .0;
            localF.resize(localF.size(), .0);
            fe.integrateDWxDN(&K, localK);
            localK *= -1;
        }
    };

public:
    enum Variable
    {
        Head,
        Velocity
    };

    void defineProblem()
    {
        //FemNodalFunctionScalar* head = new FemNodalFunctionScalar(msh, PolynomialOrder::Linear);
        //TransientFemProblem* model = new TransientFemProblem();
        _solution = new TimeEulerSpFEMLinearSolution<GWAssembly>(); 
    }

    ISolutionAlgorithm* getSolutionAlgorithm() const 
    {
        return _solution;
    }


};

TEST(Num, Discrete1)
{
    {
    // define problems
    GWFemTestProblem gwProblem;
    gwProblem.define();
    // create discrete systems according to mesh
    MeshLib::IMesh *msh = gwProblem.getMesh();
    DiscreteSystem dis(*msh, *new MathLib::SparseLinearEquations());
    // 
    gwProblem.setDiscreteSystem(dis);

    // start time stepping
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(gwProblem);

    timeStepping.setBeginning(.0);
    timeStepping.solve(100.);
    }


    //---------------------------------------------------------------------
    // How to use it
    //----------------------------------------------------------------------
    //geometry
    Rectangle rec(Point(0.0, 0.0, 0.0),  Point(2.0, 2.0, 0.0));
    MeshLib::IMesh *msh = MeshGenerator::generateStructuredRegularQuadMesh(2.0, 2, .0, .0, .0);


    GROUNDWATER_FLOW gw;
    gw.setMesh(*msh);
    gw.addBC(GROUNDWATER_FLOW::Head, *rec.getRight(), MathLib::FunctionConstant<double, double>(.0));
    gw.addST(GROUNDWATER_FLOW::Head, *rec.getLeft(), MathLib::FunctionConstant<double, double>(1.e-5));
    gw.setTransient(TimeUnit::Day, 0.0, 3.0, TimeStepFunctionConstant(0.1));
    gw.setLinearSolver(std::string("CG"), 1e-10, 1000);

    // for output
    PVD pvd(ALL_STEPS);
    pvd.addNodal(GROUNDWATER_FLOW::Head);
    pvd.addElemental(GROUNDWATER_FLOW::Velocity);
    TECPLOT tec(ALL_STEPS);
    tec.addNodal(GROUNDWATER_FLOW::Head);
    tec.setGeometry(Point(0, 0, 0));
    gw.addOutput(pvd);
    gw.addOutput(tec);

    // construct model (and check)
    gw.construct();

    // do analysis
    gw.doAnalysis();

    MathLib::IFunction<double, double>* head = gw.getResult(GROUNDWATER_FLOW::Head);
    //plot(head);

}

