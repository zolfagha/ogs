
#include <vector>

#include <gtest/gtest.h>

#include "Base/CodingTools.h"

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#include "GeoLib/Shape/Rectangle.h"

#include "MeshLib/Tools/MeshGenerator.h"
#include "MeshLib/Core/IMesh.h"

#include "FemLib/FemElementObjectContainer.h"

#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "NumLib/Discrete/ElementLocalAssembler.h"

#include "SolutionLib/IProblem.h"
#include "SolutionLib/FemProblem.h"
#include "SolutionLib/SingleStepLinearFEM.h"

#include "TestUtil.h"

using namespace GeoLib;
using namespace MathLib;
using namespace MeshLib;
using namespace FemLib;
using namespace NumLib;
using namespace SolutionLib;

class GWAssembler: public NumLib::ITimeODEElementAssembler
{
private:
    MathLib::IFunction<double, double*>* _matK;
    FemLib::LagrangianFeObjectContainer* _feObjects;
public:
    GWAssembler(FemLib::LagrangianFeObjectContainer &feObjects, MathLib::IFunction<double, double*> &mat) : _feObjects(&feObjects), _matK(&mat) 
    {
    };

    //protected:
    void assembly(const NumLib::TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::MatrixType &localM, MathLib::DenseLinearEquations::MatrixType &localK,  MathLib::DenseLinearEquations::VectorType &localF)
    {
        IFiniteElement* fe = _feObjects->getFeObject(e);

        //localM = .0;
        //localF.resize(localF.size(), .0);
        fe->integrateDWxDN(_matK, localK);
    }
};

class GWFemTestSystem : public NumLib::ITransientSystem
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





TEST(Solution, Fem1)
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

TEST(Solution, Fem2)
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
    timeStepping.solve(30.);


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

